"""Enumerate the two effective per-block source entropies under a Hadamard mask."""
import numpy as np
from scipy.special import xlogy

from run_companion_local_hadamard import project_prior


def entropy(probabilities):
    return -np.sum(xlogy(probabilities, probabilities), axis=-1)/np.log(2)


def block_entropies(p, alpha, hadamards):
    mask = np.isin(np.arange(1,10), hadamards)
    px, pz = p/(alpha+2), alpha*p/(alpha+2)
    qx = np.where(mask, px+pz, 2*px)
    qz = np.where(mask, 2*px, px+pz)
    vz_only = np.where(mask,px,pz)
    x_patterns = (np.arange(512)[:, None] >> np.arange(9)) & 1
    weights = np.prod(np.where(x_patterns, qx, 1-qx), axis=1)
    parity = np.bincount(x_patterns.sum(axis=1)%2, weights=weights)
    staged = entropy(weights)-entropy(parity)
    projected = float(entropy(project_prior(qz))[0])
    conditional = np.where(x_patterns,px/qx,vz_only/(1-qx))
    projected_exact_x = float(entropy(project_prior(conditional)) @ weights)
    return float(staged), projected, projected_exact_x


if __name__ == '__main__':
    print('p,alpha,mask,H(head|parity),H(projection),H(projection|X)')
    for p in (.07,.10,.12,.13):
        for mask in (tuple(range(1,10)), (1,3,5,7), (1,3,5,7,9)):
            values = ','.join(f'{v:.9f}' for v in block_entropies(p,100,mask))
            print(f'{p},100,{"-".join(map(str,mask))},{values}')
