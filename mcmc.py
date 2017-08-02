from emcee.utils import MPIPool
import pandas as pd
import numpy as np
import emcee

def run_emcee(run_name, nsteps, nwalkers, lnprob, to_vary):
    pool = MPIPool()
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
        
    # initiate sampler chain
    ndim = len(to_vary)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(run_name, to_vary), pool=pool)

    try:
        samples = pd.read_csv(run_name + '/samples.csv')
        print('Resuming {} at step {}'.format(run_name, samples.index[-1]//nwalkers))
        pos = np.array(samples.iloc[-nwalkers:, :-1])
        with open(run_name + '/samples.csv', 'a') as f:
            f.write('\n')
    except IOError:
        sp.call(['mkdir', run_name])
        sp.call(['mkdir', run_name + '/model_files'])
        print('Starting {}'.format(run_name))
        
        with open(run_name + '/samples.csv', 'w') as f:
            np.savetxt(f, (np.append([param[0] for param in to_vary], 'lnprob\n'),), 
                delimiter=',', fmt='%s')
        pos = [[param[1] + param[2]*np.random.randn() for param in to_vary] 
            for i in range(nwalkers)] 
            
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        print("Step {}".format(i))
        pos, chisum, blob = result
        with open(run_name + '/samples.csv', 'a') as f: 
            np.savetxt(f, [np.append(pos[i], chisum[i]) for i in range(nwalkers)], delimiter=',')

    pool.close()
    
