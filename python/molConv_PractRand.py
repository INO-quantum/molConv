#!/usr/local/bin/python3

import os
import sys
import subprocess
import io
from datetime import datetime
from struct import unpack

# test random number generators (RNG) with PractRand (0.94)
# created 21/9/2026 by Andi
# last change 21/9/2026 by Andi

# TODO: update summary which is out of date.
# TODO: finish select() for windows
# TODO: when PractRand finishs should see output of molConv on stderr but I do not see this. 
#       maybe molConv crashes or takes too long?

# - all 32bit RNG tests fail miserably!
# - Ranlux48 gives uint_fast64_t but fails test -tI/-tD indicating it cannot be used for generation of 64bit random numbers!
# - all 64bit RNG tests pass.
# - WELL1024 passes -tf test but fails -ti this is strange!? maybe has a longer period like this?
# - MersenneTwister64 should fail at about 512GB.
# - the two normally distributed tests fail at the moment since would need to convert it back into uniform distributed.
# - conclusion: do not use 32bit generators! use only 64bit generator even when floats are used in further calculation!

# summary
# test  output  RNG32   RNG48   RNG64   remark 
# -ti   32bit   fail    ok      ok      direct integer 32bit output
# -td   32bit   fail    ok      ok      float to integer 32bit conversion
# -tn   32bit   fail    fail    fail    normal distributed float to integer 32bit conversion
# -tI   64bit   fail    fail    ok      direct integer 64bit output   
# -tD   64bit   fail    fail    ok      double to integer 64bit conversion
# -tN   64bit   fail    fail    fail    normal distributed double to integer 64bit conversion
# none  -       -       -         -     self-test, stop with Ctrl-C

# test type
TEST_i = '-ti'
TEST_I = '-tI'
TEST_d = '-td'
TEST_D = '-tD'
TEST_n = '-tn'
TEST_N = '-tN'
TESTS  = [TEST_i, TEST_I, TEST_d, TEST_D, TEST_n, TEST_N]
test = '-tD'

# number of repetitions (>=1)
reps = 1

# choose data size (rounded to next higher powers of 2). MB,GB,TB are valid units
# MersenneTwister is slowest of RNG64 with 20' for 128GB
size = "8GB"

# select list of random number generators
# RNG64 are the default RNGs of molConv, all are good, MersenneTwister64 is the only from standard library but slower
# RNG24, RNG32 and RNG48 are not good and should not be used! they are here only for testing. 
#                        they must be enabled by compiling molConv with single-precision.
RNG24 = ['Ranlux24']
RNG32 = ['MinStd', 'Lehmer32', 'WELL1024']
RNG48 = ['Ranlux48']
RNG64 = ['MersenneTwister64', 'Lehmer64', 'Lehmer128', 'Wyhash64']
RNGS  = RNG24 + RNG32 + RNG48 + RNG64
test_rng = ['Lehmer64']

# molConv path and command line (%s=test, %s=generator, seed added later)
# molConv outputs molConv_ok on stderr before generating output on stdout
#molConv    = "./molConv_VisualStudio/x64/Release/molConv %s -g %s"
molConv    = "./molConv %s -g %s"
molConv_ok = 'STREAMING on stdout'

# PractRand path and command line (%s=stdin32/64, %s=size, seed added later)
# PractRand outputs PractRand_error when there is an error in command line
# PractRand outputs PractRand_fail when test(s) failed
#PractRand       = "../PractRand/PractRand_094/VisualC_net/Release/RNG_test %s -multithreaded -tlmax %s"
PractRand       = "./PractRand/PractRand_094/RNG_test %s -multithreaded -tlmax %s"
PractRand_error = 'aborting' 
PractRand_fail  = 'FAIL'
PractRand_warn  = ['unusual', 'suspicious', 'SUSPICIOUS']

# log file output (%s=RNG, %s=test)
log = "./log/%s_%s.log"

# seed to be used. None = hardware generated (default)
seed = [None, "{100,200,300,400}"][0]

# assignment of stdin32/64 with tests
STDIN64 = 'stdin64'
STDIN32 = 'stdin32'
stdin       = {TEST_i: STDIN32, TEST_I: STDIN64,
               TEST_d: STDIN32, TEST_D: STDIN64,
               TEST_n: STDIN32, TEST_N: STDIN64}
stdin_bytes = {STDIN32: 4          , STDIN64: 8}
stdin_fmt   = {STDIN32: '<I'       , STDIN64: '<Q'}

# file header also printed (%%s=date, %i/%i=repetition/repetitions, %s=size)
header  = '\n' + ('*'*80) + '\n'
header += '%%s\nPractRand test %i/%i generator %s with size = %s\n'

# colors (works only on unix systems)
if os.name == 'posix':
    color_red    = '\033[91m'
    color_green  = '\033[92m'
    color_yellow = '\033[93m'
    color_clear  = '\033[0m'
else:
    color_red    = ''
    color_green  = ''
    color_yellow = ''
    color_clear  = ''

# if not None enables test_binary_pipe to compare if first STREAM_DATA_OUT samples on stderr are the same as on stdout
STREAM_DATA_OUT = None

def test_binary_pipe(psrc):
    # on Windows we need to set pipe in binary mode, otherwise newlines are inserted into the data which causes all tests to fail!
    # in C++ code the option STREAM_DATA_OUT can be enabled and it outputs first STREAM_DATA_OUT samples in parallel on stderr and stdout.
    # here we print first STREAM_DATA_OUT samples from molConv and compare what molConv prints on stderr.
    # when they are the same then pipe is working correctly in binary mode, otherwise not. test several times to be sure!
    # notes: 
    # - in 2023/2024 (Windows 8 or 10) I managed to get binary pipe with "cmd /c runtests.bat" 
    #   with the .bat file launching molConv and piping stdout into PractRand.
    #   this does not work anymore in 2026 on Windows 10!
    # - using in python msvcrt.setmode does not work either. return value indicates that previous mode is binary but which is wrong!
    # - the new solution is to use _setmode in molConv C++ code.

    # try to set binary mode (no error but does not work)
    if os.name == 'nt':
        import msvcrt
        mode = msvcrt.setmode(psrc.stdout.fileno(), os.O_BINARY)
        print('stdout set mode to binary, old mode is binary:', mode == os.O_BINARY)

    # read STREAM_DATA_OUT samples from stdout
    bytes_to_read = stdin_bytes[test_stdin]
    fmt           = stdin_fmt  [test_stdin]
    count = 0
    data = []
    while psrc.poll() is None and count < STREAM_DATA_OUT:
        line = psrc.stdout.read(bytes_to_read)
        if len(line) == 0:
            break
        value = unpack(fmt, line)
        if bytes_to_read == 8:
            print('molConv >>> 0x%016x' % (value))
            data.append('0x%016x' % (value))
        else:
            print('molConv >>> 0x%08x' % (value))
            data.append('0x%08x' % (value))        
        count += 1
    psrc.stdout.close() # causes molConv to close and return
    print()
    
    # read STREAM_DATA_OUT lines from stderr
    # this reads empty lines when psrc is closed
    count = 0
    fail_count = 0
    while count < STREAM_DATA_OUT:
        line = psrc.stderr.readline()
        if len(line) == 0:
            break
        line = line.strip().decode('utf-8')
        if len(line) == (bytes_to_read*2+2) and line.startswith('0x'):
            if line != data[count]:
                print('molConv   > %s != %s error!' % (line, data[count]))
                fail_count += 1
            else:
                print('molConv   > %s ok' % line)
            count += 1
        elif len(line) > 0:
            print('molConv   >', line)
    if count == STREAM_DATA_OUT:
        if fail_count == 0:
            print('%i/%i random numbers received ok\n' % (count, STREAM_DATA_OUT))
        else:
            print('%i/%i random numbers received with errors!\n' % (count, STREAM_DATA_OUT))
    elif count == 0:
        if test == '-tI' or test == '-ti':
            print('%i/%i random numbers received! is STREAM_DATA_OUT enabled in molConv C++ code?\n' % (count, STREAM_DATA_OUT))
        else:
            print("%i/%i random numbers received! test must be '-tI' or '-ti'! actual test = '%s'\n" % (count, STREAM_DATA_OUT, test))
    else:
        print('%i/%i random numbers received!? %i failed.\n' % (count, STREAM_DATA_OUT, fail_count))

if os.name == 'nt':
    # on Windows select does not work with stdin, so we have to do it manually
    # TODO: this is not finished!
    import ctypes
    def select(rlist, wlist, elist, timeout=None):
        # get windows handles
        _hr = [r.fileno() for r in rlist]
        _hw = [w.fileno() for w in wlist]
        _he = [e.fileno() for e in elist]
        # wait for handles
        _h = _hr + _hw + _he
        arrtype = ctypes.c_long * len(_h)
        _ph = arrtype(_h) 
        result = ctypes.windll.kernel32.WaitForMultipleObjects(len(_h), _ph, False, int(timeout*1000) if timeout is not None else -1)
        if result == 0x102:
            print('timeout')
            return [[],[],[]]
        elif result == 0xffffffff:
            print('error')
            return [[],[],[]]
        else:
            h = _h[result]
            if   result < len(_hr):          return [[h],[],[]]
            elif result < len(_hr)+len(_hw): return [[],[h],[]]
            else:                            return [[],[],[h]]
        return result
else:
    from select import select

if __name__ == '__main__':

    try:
        test_stdin = stdin[test]
    except KeyError:
        print("unknown test '%s'\n" %  test)
        exit()
        
    if os == 'nt':
        # required on Windows, otherwise process cannot be killed and runs forever!
        startupinfo = subprocess.CREATE_NEW_PROCESS_GROUP
    else:
        startupinfo = None

    fail_count = 0
    total_count = 0
    for rng in test_rng:
        for rep in range(reps):
            t_start = datetime.now()
            ok = True
            info = header % (rep+1, reps, rng, size)
            _info = info % (str(t_start))
            
            if rng not in RNGS:
                print("unknown RNG '%s'\n" % rng)
                exit()
                        
            src = molConv % (test, rng)
            dst   = PractRand % (test_stdin, size)
            if seed is not None:
                src += ' -s ' + seed
                        
            filename = log % (rng, test)
            
            with open(filename, "ab+") as f:
            
                f.write(_info.encode('utf-8') + b'\n')
                print(_info)
                print('log file  :', filename)
                print('molConv   :', src)
                print('PractRand :', dst)
                print()
    
                # start molConv and wait for seed values written on stderr
                
                psrc = subprocess.Popen(src.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, startupinfo=startupinfo)

                seed_act = None
                while psrc.poll() is None:
                    line = psrc.stderr.readline()
                    if len(line) == 0:
                        break
                    else:
                        f.write(line)
                        line = line.strip().decode('utf-8')
                        start = line.find('{')
                        if start >= 0:
                            end = line.find('}', start)
                            seed_act = line[start:start + end]
                            print('molConv   >', line[:start] + color_green + seed_act + color_clear)
                        else:
                            print('molConv   >', line)
                        if molConv_ok in line:
                            break
                print()
                
                if psrc.poll() is not None:
                    print('error in molConv! please check output!\n')
                    exit()

                # check if we have got seed values so we can pass them to PractRand and save to log file
                if seed_act is not None:
                    dst += ' -seed ' + seed_act
                    #f.write(('seed = ' + seed_act + '\n\n').encode('utf-8'))

                if STREAM_DATA_OUT is not None:
                    test_binary_pipe(psrc)
                    exit()
                    
                # start PractRand with stdin from stdout of molConv
                pdst = subprocess.Popen(dst.split(' '), stdout=subprocess.PIPE, stdin=psrc.stdout, startupinfo=startupinfo)

                # wait for output of PractRand and molConv until they finish
                #if os.name == 'nt':
                if True:
                    # select does not work on Windows. fallback to blocking calls on individual processes
                    # wait for PractRand to close
                    while True:
                        line = pdst.stdout.readline()
                        if len(line) == 0: 
                            print('PractRand : stdout closed!')
                            psrc.stdout.close() # this tells molConv to close
                            break
                        else:
                            f.write(line)
                            line = line.strip().decode('utf-8')
                            if PractRand_fail in line or PractRand_error in line:
                                ok = False
                                color = [color_red, color_clear]
                            else: 
                                color = ['','']
                                for w in PractRand_warn:
                                    if w in line:
                                        color = [color_yellow, color_clear]
                                        break
                            print(color[0] + 'PractRand >', line + color[1])
                    # wait for molConv to close
                    # TODO: this blocks with other tests than -tI, -ti which output data on stderr
                    if test == TEST_I or test == TEST_i:
                        while True:
                            line = psrc.stderr.readline()
                            if len(line) == 0: 
                                print('molConv   : stderr closed!')
                                psrc.stdout.close() # this should close PractRand when still running
                                break
                            else:
                                line = line.strip().decode('utf-8')
                                print('molConv   >', line)
                else:
                    rd = [pdst.stdout, psrc.stderr]
                    while len(rd) > 0:
                        result = select(rd, [], [], 1)
                        if len(rd) == 1 and len(result[0]) == 0:
                            print('timeout with single process running!')
                            break
                        if pdst.stdout in result[0]:
                            line = pdst.stdout.readline()
                            if len(line) == 0: 
                                print('PractRand : stdout closed!')
                                rd = [r for r in rd if r != pdst.stdout]
                                psrc.stdout.close() # this tells molConv to close
                            else:
                                f.write(line)
                                line = line.strip().decode('utf-8')
                                if PractRand_fail in line or PractRand_error in line:
                                    ok = False
                                    color = [color_red, color_clear]
                                else: 
                                    color = ['','']
                                    for w in PractRand_warn:
                                        if w in line:
                                            color = [color_yellow, color_clear]
                                            break
                                print(color[0] + 'PractRand >', line + color[1])
                        if psrc.stderr in result[0]:
                            line = psrc.stderr.readline()
                            if len(line) == 0: 
                                print('molConv   : stderr closed!')
                                rd = [r for r in rd if r != psrc.stderr]
                                psrc.stdout.close() # this should close PractRand when still running
                            else:
                                line = line.strip().decode('utf-8')
                                print('molConv   >', line)
                   
                # ensure molConv and PractRand are closed, otherwise kill them.
                # TODO: on Windows processes must be killed?
                if psrc.poll() is None:
                    print('molConv   : killed!')
                    psrc.kill()
                else:
                    print('molConv   : closed. return code', psrc.returncode)
                if pdst.poll() is None:
                    print('PractRand : killed!')
                    pdst.kill()
                else:
                    print('PractRand : closed. return code', pdst.returncode)

                _info = (info % (str(datetime.now()))) + 'done in ' + ('%.1f'%(datetime.now()-t_start).seconds) + ' seconds, result = '
                result = 'SUCCESS' if ok else 'FAIL!'
                
                f.write((_info+result+'\n').encode('utf-8'))
                color = color_green if ok else color_red
                print(_info + color + result + color_clear)

            if not ok: fail_count += 1
            total_count += 1

    if fail_count == 0:
        print(color_green + ('\ndone %i tests ok\n' % (total_count)) + color_clear)
    else:
        print(color_red + '\ndone %i/%i tests failed!\n' % (fail_count, total_count) + color_clear)

