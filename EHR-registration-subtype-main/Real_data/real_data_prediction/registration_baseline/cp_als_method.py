# extracted from pyTen package on github
import numpy as np
import torch
from silrtc import permute

def khatrirao(u):
    """
    Calculate The Khatrirao Product, param u: a list of 2-D arrays
    """
    r = u[0].shape[1]
    k = []
    n = len(u)
    for j in range(r):
        temp = 1
        for i in range(n):
            temp1 = np.outer(temp, u[i][:,j])
            temp = temp1.reshape([1, temp1.size])
        k = np.append(k, temp)
    k = (k.reshape([r, int(len(k)/r)])).transpose()
    return k

class Ktensor(object):
    """
    Tensor stored in decomposed form as a Kruskal operator (CP decomposition).
    """
    def __init__(self, lmbda=None, us=None):
        """
        Constructor for Ktensor (CP Tensor) object with the weights and latent matrices.
        ----------
        :type self: object
        :param lmbda : array_like of floats, optional
           Weights for each dimension of the Kruskal operator.
           ``len(lambda)`` must be equal to ``U[i].shape[1]``
        :param us : list of ndarrays
           Factor matrices from which the Tensor representation
           is created. All factor matrices ``U[i]`` must have the
           same number of columns, but can have different
           number of rows.
        ----------
        """
        if us is None:
            raise ValueError("Ktensor: first argument cannot be empty.")
        else:
            self.Us = np.array(us)
        self.shape = tuple(Ui.shape[0] for Ui in us)
        self.ndim = len(self.Us)
        self.rank = self.Us[0].shape[1]
        if lmbda is None:
            self.lmbda = np.ones(len(self.rank))
        else:
            self.lmbda = np.array(lmbda)
        if not all(np.array([Ui.shape[1] for Ui in us]) == self.rank):
            raise ValueError('Ktensor: dimension mismatch of factor matrices')

    def norm(self):
        """
        Efficient computation of the Frobenius norm for ktensors
        Returns: None
        -------
        norm : float
               Frobenius norm of the Ktensor
        """
        coefmatrix = np.dot(self.Us[0].T, self.Us[0])
        for i in range(1, self.ndim):
            coefmatrix = coefmatrix * np.dot(self.Us[i].T, self.Us[i])
        coefmatrix = np.dot(np.dot(self.lmbda.T, coefmatrix), self.lmbda)
        return np.sqrt(coefmatrix.sum())

    def tondarray(self):
        """
        Converts a Ktensor into a dense multidimensional ndarray
        Returns: None
        -------
        arr : np.ndarray
            Fully computed multidimensional array whose shape matches
            the original Ktensor.
        """
        a = np.dot(self.lmbda.T, khatrirao(self.Us).T)
        return a.reshape(self.shape)

def unfold(x, n=None):
        """Return the mode-n unfold of the Tensor."""
        if n is None:
            raise ValueError('Tensor/UNFOLD: unfold mode n (int) needs to be specified.')
        N = x.ndim
        temp1 = [n]
        temp2 = range(n)
        temp3 = range(n+1, N)
        temp1[len(temp1):len(temp1)] = temp2
        temp1[len(temp1):len(temp1)] = temp3
        xn = permute(x, temp1)
        xn = xn.reshape([xn.shape[0], np.prod(xn.shape)//xn.shape[0]])
        return xn

def nvecs(x, n=None, r=None):
        """Return first r eigenvectors of the mode-n unfolding matrix"""
        if n is None:
            raise ValueError('Tensor/NVECS: unfold mode n (int) needs to be specified.')
        if r is None:
            raise ValueError('Tensor/NVECS: the number of eigenvectors r needs to be specified.')
        xn = unfold(x, n)
        [eigen_value, eigen_vector] = np.linalg.eig(xn.dot(xn.transpose()))
        return eigen_vector[:, range(r)]

def cp_als(y, r=10, omega=None, tol=1e-3, maxiter=10, init='random', printitn=10):
    """ CP_ALS Compute a CP decomposition of a Tensor (and recover it).
    ---------
     :param  'y' - Tensor with Missing data
     :param  'r' - Rank of the tensor
     :param 'omega' - Missing data Index Tensor
     :param 'tol' - Tolerance on difference in fit
     :param 'maxiters' - Maximum number of iterations
     :param 'init' - Initial guess ['random'|'nvecs'|'eigs']
     :param 'printitn' - Print fit every n iterations; 0 for no printing
    ---------
     :return
        'P' - Decompose result.(kensor)
        'X' - Recovered Tensor.
    ---------
    """
    X = np.copy(y)
    if omega is None:
        omega = X * 0 + 1

    # Extract number of dimensions and norm of X.
    N = X.ndim
    X = np.nan_to_num(X)
    normX = torch.norm(torch.Tensor(X)).item()
    dimorder = list(range(N)) # 'dimorder' - Order to loop through dimensions {0:(ndims(A)-1)}

    # Define convergence tolerance & maximum iteration
    fitchangetol = tol
    maxiters = maxiter

    # Recover or just decomposition
    recover = 0
    if 0 in omega:
        recover = 1

    # Set up and error checking on initial guess for U.
    if type(init) == list:
        Uinit = init[:]
        if len(Uinit) != N:
            raise IndexError('OPTS.init does not have %d lists', N)
        for n in dimorder[1:]:
            if Uinit[n].shape != (X.shape[n], r):
                raise IndexError('OPTS.init{%d} is the wrong size', n)
    else:
        # Observe that we don't need to calculate an initial guess for the
        # first index in dimorder because that will be solved for in the first
        # inner iteration.
        if init == 'random':
            Uinit = [[] for _ in range(N)]
            #Uinit[0] = []
            for n in dimorder[1:]:
                Uinit[n] = np.random.random([X.shape[n], r])
        elif init == 'nvecs' or init == 'eigs':
            Uinit = range(N)
            Uinit[0] = []
            for n in dimorder[1:]:
                Uinit[n] = nvecs(X, n, r) # first r leading eigenvecters
        else:
            raise TypeError('The selected initialization method is not supported')

    # Set up for iterations - initializing U and the fit.
    U = Uinit[:]
    fit = 0

    if printitn > 0:
        print('\nCP_ALS:\n')

    # Save hadamard product of each U[n].T*U[n]
    UtU = np.zeros([N, r, r])
    for n in range(N):
        if len(U[n]):
            UtU[n,:,:] = np.dot(U[n].T, U[n])

    for iter in range(1, maxiters+1):
        fitold = fit
        oldX = X * 1.0

        # Iterate over all N modes of the Tensor
        for n in range(N):
            # Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
            temp1 = [n]
            temp2 = list(range(n))
            temp3 = list(range(n+1, N))
            temp2.reverse()
            temp3.reverse()
            temp1[len(temp1):len(temp1)] = temp3
            temp1[len(temp1):len(temp1)] = temp2
            Xn = permute(X, temp1)
            Xn = Xn.reshape([Xn.shape[0], int(Xn.size/Xn.shape[0])])
            tempU = U[:]
            tempU.pop(n)
            tempU.reverse()
            Unew = Xn.dot(khatrirao(tempU))

            # Compute the matrix of coefficients for linear system
            temp = list(range(n))
            temp[len(temp):len(temp)] = list(range(n+1, N))
            y = np.prod(UtU[temp,:,:], axis=0)
            Unew = Unew.dot(np.linalg.inv(y))

            # Normalize each vector to prevent singularities in coefmatrix
            if iter == 1:
                lamb = np.sqrt(np.sum(np.square(Unew), 0)) # 2-norm
            else:
                lamb = np.max(Unew, 0)
                lamb = np.max([lamb, np.ones(r)], 0) # max-norm

            lamb = [x * 1.0 for x in lamb]
            Unew = Unew / np.array(lamb)
            U[n] = Unew
            UtU[n, :, :] = np.dot(U[n].T, U[n])

        # Reconstructed fitted Ktensor
        P = Ktensor(lamb, U)
        if recover == 0:
            if normX == 0:
                fit = P.norm() ** 2 - 2 * np.sum(X * P.tondarray())
            else:
                normresidual = np.sqrt(abs(normX ** 2 + P.norm() ** 2 - 2 * np.sum(X * P.tondarray())))
                fit = 1 - (normresidual / normX) # fraction explained by model
                fitchange = abs(fitold - fit)
        else:
            temp = P.tondarray()
            X = X * (1 - omega) + temp * omega
            fitchange = np.linalg.norm(X - oldX)

        # Check for convergence
        if (iter > 1) and (fitchange < fitchangetol):
            flag = 0
        else:
            flag = 1

        if ((printitn!=0) and (iter % printitn==0)) or ((printitn>0) and (flag==0)):
            if recover == 0:
                print('CP_ALS: iterations={0}, f={1}, f-delta={2}'.format(iter, fit, fitchange))
            else:
                print('CP_ALS: iterations={0}, f-delta={1}'.format(iter, fitchange))

        # Check for convergence
        if flag == 0:
            break
    return P, X
