import numpy as np
import math
import torch

def permute(x, order=None):
    """ Returns a Tensor permuted by the order specified."""
    if order.__class__ == list or order.__class__ == tuple:
        order = np.array(order)
    if x.ndim != len(order):
        raise ValueError("Permute: Invalid permutation order.")
    if not (sorted(order) == np.arange(x.ndim)).all():
        raise ValueError("Permute: Invalid permutation order.")
    newdata = np.copy(x)
    newdata = newdata.transpose(order)
    return newdata

def prod(arg):
    """ returns the product of elements in arg.
    arg can be list, tuple, set, and array with numerical values. """
    ret = 1
    for i in range(0, len(arg)):
        ret = ret * arg[i]
    return ret

def ipermute(x, order=None):
    """ Returns a Tensor permuted by the inverse of the order specified """
    if order.__class__ == np.array or order.__class__ == tuple:
        order = list(order)
    else:
        if order.__class__ != list:
            raise ValueError('Ipermute: permutation order must be a list.')

    if x.ndim != len(order):
        raise ValueError("Ipermute: invalid permutation order.")
    if not (sorted(order) == np.arange(x.ndim)).all():
        raise ValueError("Ipermute: invalid permutation order.")
    iorder = [order.index(idx) for idx in range(0, len(order))]
    return permute(x, iorder)

class Tenmat(object):
    """
    Store a Matricization of a Tensor object.
    """
    def __init__(self, x=None, rdim=None, cdim=None, tsize=None):
        """
        Create a Tenmat object from a given Tensor X
         ----------
        :param x: dense Tensor object.
        :param rdim: an one-dim array representing the arranged dimension index for the matrix column
        :param cdim: an one-dim array representing the arranged dimension index for the matrix row
        :param tsize: a tuple denoting the size of the original tensor
        :return: constructed Matricization of a Tensor object.
        ----------
        """
        # convert a Tensor to a matrix
        if rdim.__class__ == list or rdim.__class__ == int:
            rdim = np.array(rdim) - 1

        self.shape = x.shape

        if cdim is None:
            cdim = np.array([y for y in range(0, x.ndim) if y not in np.zeros(x.ndim - 1) + rdim])
        elif cdim.__class__ == list or cdim.__class__ == int:
            cdim = np.array(cdim) - 1
        else:
            raise ValueError("Tenmat: incorrect specification of dimensions.")

        if not (np.arange(x.ndim) == sorted(np.append(rdim, cdim))).all():
            raise ValueError("Tenmat: second argument must be a list or an integer.")

        self.rowIndices = rdim
        self.colIndices = cdim
        x = permute(x, np.append(rdim, cdim))

        if type(rdim) != np.ndarray:
            row = prod([self.shape[y] for y in [rdim]])
        else:
            row = prod([self.shape[y] for y in rdim])
        if type(cdim) != np.ndarray:
            col = prod([self.shape[y] for y in [cdim]])
        else:
            col = prod([self.shape[y] for y in cdim])

        self.data = x.reshape([row, col], order='F')

    def totensor(self):
        # returns a Tensor object based on a Tenmat
        order = np.append(self.rowIndices, self.colIndices)
        data = self.data.reshape([self.shape[idx] for idx in order], order='F')
        t_data = ipermute(data, list(order))
        return t_data

def pro_to_trace_norm(z, tau):
    m = z.shape[0]
    n = z.shape[1]
    if 2 * m < n:
        [U, Sigma2, V] = np.linalg.svd(np.dot(z, z.T))
        S = np.sqrt(Sigma2)
        tol = np.max(z.shape) * (2 ** int(math.log(max(S), 2))) * 2.2204 * 1E-16
        k = np.sum(S > max(tol, tau))
        mid = [max(S[i] - tau, 0) * 1.0 / S[i] for i in range(k)]
        X = np.dot(np.dot(U[:, 0:k], np.dot(np.diag(mid), U[:, 0:k].T)), z)
        return X, k, Sigma2

    if m > 2 * n:
        z = z.T
        [U, Sigma2, V] = np.linalg.svd(np.dot(z, z.T))
        S = np.sqrt(Sigma2)
        tol = np.max(z.shape) * (2 ** int(math.log(max(S), 2))) * 2.2204 * 1E-16
        k = np.sum(S > max(tol, tau))
        mid = [max(S[i] - tau, 0) * 1.0 / S[i] for i in range(k)]
        X = np.dot(np.dot(U[:, 0:k], np.dot(np.diag(mid), U[:, 0:k].T)), z)
        return X.T, k, Sigma2

    [U, S, V] = np.linalg.svd(z)
    Sigma2 = S ** 2
    k = sum(S > tau)
    X = np.dot(U[:, 0:k], np.dot(np.diag(S[0:k] - tau), V[0:k, :]))
    return X, n, Sigma2

def silrtc_model(x, omega=None, alpha=None, gamma=None, max_iter=10, epsilon=1e-3, printitn=10):
    """
    Simple Low Rank Tensor Completion (SiLRTC).
    Reference: "Tensor Completion for Estimating Missing Values in Visual Data", PAMI, 2012.
    """
    # x and omega must be numpy ndarray!
    T = np.copy(x)
    N = x.ndim
    if printitn == 0:
        printitn = max_iter
    if omega is None:
        omega = x * 0 + 1
    if alpha is None:
        alpha = np.ones([N])
        alpha = alpha / sum(alpha)
    if gamma is None:
        gamma = 0.1 * np.ones([N])

    # initialization
    x = np.nan_to_num(x)
    normX = torch.norm(torch.Tensor(x)).item()
    print(normX)
    errList = np.zeros([max_iter, 1])

    #M = np.arange(N)
    gammasum = sum(gamma)
    tau = alpha / gamma
    for k in range(max_iter):
        print('Epoch: '+str(k))

        if ((k + 1) % printitn == 0) and (k != 0) and (printitn != max_iter):
            print('SiLRTC: iterations = {0} difference = {1}\n'.format(k, errList[k - 1]))

        Xsum = 0
        for i in range(N):
            temp = Tenmat(x, i+1)
            temp1, tempn, tempSigma2 = pro_to_trace_norm(temp.data, tau[i])
            temp.data = temp1
            M = temp.totensor()
            Xsum = Xsum + gamma[i] * M

        Xlast = np.copy(x)
        x = Xsum / gammasum
        x = x * omega + np.nan_to_num(T) * (1 - omega)
        #print(np.sum(np.isnan(x)))
        diff = x - Xlast
        errList[k] = np.linalg.norm(diff) / normX
        if errList[k] < epsilon:
            errList = errList[0:(k + 1)]
            break

    print('SiLRTC ends: total iterations = {0} difference = {1}\n\n'.format(k + 1, errList[k]))
    return x
