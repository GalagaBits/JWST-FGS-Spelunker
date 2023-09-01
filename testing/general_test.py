"""General testing suite for spelunker

This script is intended to be run before pushing any release of the code. Passing means the code works only, not that the code is 
producing results that make sense.

"""

import os
import shutil

import spelunker

def run_tests():

    results = {}
    results['download'] = {}
    results['download']['description'] = 'Downloading & storage of data (spk.load)'
    results['download']['pass'] = False
    results['download']['function'] = test_download
    results['2dgaussian'] = {}
    results['2dgaussian']['pass'] = False
    results['2dgaussian']['function'] = test_2dgaussian
    results['2dgaussian']['description'] = '2D gaussian fits (spk.gauss2d_fit)'

    spk = 'init'
    for key in list( results.keys() ):

        results[key]['pass'], spk, results[key]['errorname'], results[key]['error'] = results[key]['function'](spk)

        if not results[key]['pass']:

            break

    test_summary(results)

def test_download(spk):

    try: 

        spk = spelunker.load(pid=1541, obs_num='1', visit='1', save=True)

        return True, spk, '', ''

    except Exception as error:

        return False, None, type(error).__name__, error

def test_2dgaussian(spk):

    # Fits to only first 100 frames:
    spk.fg_array = spk.fg_array[:100]

    try:

        spk.gauss2d_fit()

        return True, spk, '', ''

    except Exception as error:

        return False, None, type(error).__name__, error

def test_summary(results):

    all_good = True
    print('\n >> Results:\n')
    for key in list( results.keys() ):

        if results[key]['pass']:

            print('\t - '+results[key]['description']+'........\033[32mPASS\033[0m')

        else:

            print('\t - '+results[key]['description']+'........\033[31mFAIL\033[0m')
            all_good = False
            
    if not all_good:

        print('\n\t - Final assessment: some tests \033[31mfailed\033[0m. Here is a summary of those, along with errors raised:\n')

        counter = 1
        for key in list( results.keys() ):

            if not results[key]['pass']:

                try:

                    print('\t\t '+str(counter)+') '+results[key]['description']+' \033[31mfailed\033[0m with a',results[key]['errorname'],':\n')
                    print(results[key]['error'])
                    print('\n')

                except:

                    print('\t\t '+str(counter)+') '+results[key]['description']+' test not run (possibily because of dependance with prior test that failed)')

                counter += 1
        print('\n >> Some spelunker tests \033[31mfailed\033[0m. Please check errors above.\n')

    else:

        print('\n\t - Final assessment: all tests \033[32msuccessful\033[0m!')

        print('\n >> spelunker is \033[32mready to go\033[0m!\n')

print('\n >> Testing spelunker version ',spelunker.__version__)

if os.path.exists('spelunker_outputs'):

    shutil.rmtree('spelunker_outputs')

run_tests()
