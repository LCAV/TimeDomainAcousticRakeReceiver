
import numpy as np
from scipy.signal import resample,fftconvolve
from scipy.linalg import toeplitz

import pyroomacoustics as pra


'''
We create a new Beamformer class for Rake Perceptually motivated beamformers
'''
class RakePerceptual_TD(pra.Beamformer):

    def computeWeights(self, sources, interferers, R_n, epsilon=1e-2):

        dist_mat = pra.distance(self.R, sources)
        s_time = dist_mat / pra.c
        s_dmp = 1./(4*np.pi*dist_mat)

        dist_mat = pra.distance(self.R, interferers)
        i_time = dist_mat / pra.c
        i_dmp = 1./(4*np.pi*dist_mat)

        # compute offset needed for decay of sinc by epsilon
        offset = np.maximum(s_dmp.max(), i_dmp.max())/(np.pi*self.Fs*epsilon)
        t_min = np.minimum(s_time.min(), i_time.min())
        t_max = np.maximum(s_time.max(), i_time.max())
        
        kappa = int(t_max*float(self.Fs))
        kappa = int(0.03*self.Fs)

        # adjust timing
        s_time -= t_min - offset
        i_time -= t_min - offset
        Lh = int((t_max - t_min + 2*offset)*float(self.Fs))

        # the channel matrix
        K = sources.shape[1]
        Lg = self.Lg
        off = (Lg - Lh)/2
        L = self.Lg + Lh - 1

        H = np.zeros((Lg*self.M, 2*L))

        for r in np.arange(self.M):

            # build interferer RIR matrix
            hx = pra.lowPassDirac(s_time[r,:,np.newaxis], s_dmp[r,:,np.newaxis], self.Fs, Lh).sum(axis=0)
            #row = np.pad(hx, ((0,L-len(hx))), mode='constant')
            #col = np.pad(hx[:1], ((0, Lg-1)), mode='constant')
            #H[r*Lg:(r+1)*Lg,:L] = toeplitz(col, row)
            H[r*Lg:(r+1)*Lg,:L] = pra.convmtx(hx, Lg).T

            # build interferer RIR matrix
            hq = pra.lowPassDirac(i_time[r,:,np.newaxis], i_dmp[r,:,np.newaxis], self.Fs, Lh).sum(axis=0)
            #row = np.pad(hq, ((0,L-len(hq))), mode='constant')
            #col = np.pad(hq[:1], ((0, Lg-1)), mode='constant')
            #H[r*Lg:(r+1)*Lg,L:] = toeplitz(col, row)
            H[r*Lg:(r+1)*Lg,L:] = pra.convmtx(hq, Lg).T

        # the constraint vector
        c = np.zeros(L)
        n0 = kappa
        n1 = np.minimum(n0 + np.ceil(0.03*self.Fs), L)
        n2 = np.minimum(n0 + np.ceil(0.05*self.Fs), L)
        exponent = -1./(0.02*self.Fs)*np.log(1e-2)
        c[n0:n1] = np.ones(n1-n0)
        if (n2 > n1):
            c[n1:n2] = np.exp(-exponent*np.arange(n2-n1))
        c = c[:,np.newaxis]

        A = H[:,:n0+1]
        b = np.zeros((n0+1,1))
        b[-1,0] = 1

        h = H[:,n0,np.newaxis]

        # We first assume the sample are uncorrelated
        K_nq = np.dot(H[:,L:], H[:,L:].T) + R_n

        # causal response construction
        B = np.dot(np.linalg.inv(K_nq), A)
        C = np.linalg.inv(np.dot(A.T, B))
        g_val = np.dot(B, np.dot(C, b))

        '''
        from cvxopt import matrix, solvers

        # rename everything in QP terminology
        Q = 2.*matrix(K_nq)
        p = matrix(np.zeros(self.M*Lg))
        G = matrix(np.concatenate((H[:,:n0], H[:,n0+1:L]), axis=1).T)
        c1 = np.concatenate((c[:n0,:], c[n0+1:,:]), axis=0)
        h = matrix(c1[:,0])
        A = matrix(H[:,n0], (1,self.M*Lg))
        b = matrix(1.0)
        
        # solve the QP
        sol = solvers.qp(Q, p, G, h, A, b)
        g_val = np.array(sol['x'])
        '''

        # reshape and store
        self.filters = g_val.reshape((self.M, self.Lg))

        '''
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(np.arange(L)/float(self.Fs), np.dot(H[:,:L].T, g_val))
        plt.plot(np.arange(L)/float(self.Fs), np.dot(H[:,L:].T, g_val))
        plt.legend(('Channel of desired source','Channel of interferer'))
        '''

        A = np.dot(g_val.T, H[:,:L])
        num = np.dot(A, A.T)
        denom =  np.dot(np.dot(g_val.T, K_nq), g_val)
        return num/denom



'''
We create a new Beamformer class for Rake MaxSINR in time-domain
'''
class RakeMaxSINR_TD(pra.Beamformer):

    def computeWeights(self, sources, interferers, R_n, epsilon=1e-2):

        dist_mat = pra.distance(self.R, sources)
        s_time = dist_mat / pra.c
        s_dmp = 1./(4*np.pi*dist_mat)

        dist_mat = pra.distance(self.R, interferers)
        i_time = dist_mat / pra.c
        i_dmp = 1./(4*np.pi*dist_mat)

        # compute offset needed for decay of sinc by epsilon
        offset = np.maximum(s_dmp.max(), i_dmp.max())/(np.pi*self.Fs*epsilon)
        t_min = np.minimum(s_time.min(), i_time.min())
        t_max = np.maximum(s_time.max(), i_time.max())

        # adjust timing
        s_time -= t_min - offset
        i_time -= t_min - offset
        Lh = int((t_max - t_min + 2*offset)*float(self.Fs))

        # the channel matrix
        K = sources.shape[1]
        Lg = self.Lg
        off = (Lg - Lh)/2
        L = self.Lg + Lh - 1

        H = np.zeros((Lg*self.M, 2*L))

        for r in np.arange(self.M):

            # build interferer RIR matrix
            hx = pra.lowPassDirac(s_time[r,:,np.newaxis], s_dmp[r,:,np.newaxis], self.Fs, Lh).sum(axis=0)
            H[r*Lg:(r+1)*Lg,:L] = pra.convmtx(hx, Lg).T

            # build interferer RIR matrix
            hq = pra.lowPassDirac(i_time[r,:,np.newaxis], i_dmp[r,:,np.newaxis], self.Fs, Lh).sum(axis=0)
            H[r*Lg:(r+1)*Lg,L:] = pra.convmtx(hq, Lg).T

        # the constraint vector
        c = np.zeros(L)
        kappa = int(0.03*self.Fs)
        c[kappa] = 1.

        # We first assume the sample are uncorrelated
        K_nq = np.dot(H[:,L:], H[:,L:].T) + R_n

        '''
        # Compute the TD filters
        K_nq_inv = np.linalg.inv(K_nq)
        C = np.dot(K_nq_inv, H[:,:L])
        B = np.linalg.inv(np.dot(H[:,:L].T, C))
        g_val = np.dot(C, np.dot(B, c))
        self.filters = g_val.reshape((self.M,Lg))
        '''

        #'''
        # faster (?)
        A = pra.levinson(K_nq[:,0], H[:,:L])
        B = H[:,:L].T.dot(A)
        x = np.linalg.solve(B, c)
        g_val = np.dot(A, x)
        self.filters = g_val.reshape((self.M,Lg))
        #'''

        '''
        # Compute TD filters using generalized Rayleigh coefficient maximization
        C_inv = np.linalg.inv(np.linalg.cholesky(K_nq))
        B = np.dot(H[:,:L].T, C_inv)
        l, v = np.linalg.eig( np.dot(B.T,B) ) 
        print v[:,0].sum()
        g_val = np.dot(C_inv.T, np.real(v[:,0]))
        self.filters = g_val.reshape((self.M, Lg))
        '''

        '''
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(np.arange(L)/float(self.Fs), np.dot(H[:,:L].T, g_val))
        plt.plot(np.arange(L)/float(self.Fs), np.dot(H[:,L:].T, g_val))
        plt.legend(('Channel of desired source','Channel of interferer'))
        '''

        A = np.dot(g_val.T, H[:,:L])
        num = np.dot(A, A.T)
        denom =  np.dot(np.dot(g_val.T, K_nq), g_val)
        return num/denom


'''
We create a new Beamformer class for Rake One-Forcing in time-domain
'''
class RakeOF_TD(pra.Beamformer):

    def computeWeights(self, sources, interferers, R_n, epsilon=1e-2):

        dist_mat = pra.distance(self.R, sources)
        s_time = dist_mat / pra.c
        s_dmp = 1./(4*np.pi*dist_mat)

        dist_mat = pra.distance(self.R, interferers)
        i_time = dist_mat / pra.c
        i_dmp = 1./(4*np.pi*dist_mat)

        # compute offset needed for decay of sinc by epsilon
        offset = np.maximum(s_dmp.max(), i_dmp.max())/(np.pi*self.Fs*epsilon)
        t_min = np.minimum(s_time.min(), i_time.min()) - offset
        t_max = np.maximum(s_time.max(), i_time.max()) + offset

        # adjust timing
        s_time -= t_min
        i_time -= t_min
        Lh = int((t_max - t_min)*float(self.Fs))

        # the channel matrix
        K = sources.shape[1]
        Lg = self.Lg
        off = (Lg - Lh)/2
        L = self.Lg + Lh - 1

        H = np.zeros((Lg*self.M, 2*L))
        As = np.zeros((Lg*self.M, K))

        for r in np.arange(self.M):

            # build constraint matrix
            hs = pra.lowPassDirac(s_time[r,:,np.newaxis], s_dmp[r,:,np.newaxis], self.Fs, Lh)[:,::-1]
            As[r*Lg+off:r*Lg+Lh+off,:] = hs.T

            # build interferer RIR matrix
            hx = pra.lowPassDirac(s_time[r,:,np.newaxis], s_dmp[r,:,np.newaxis], self.Fs, Lh).sum(axis=0)
            #row = np.pad(hx, ((0,L-len(hx))), mode='constant')
            #col = np.pad(hx[:1], ((0, Lg-1)), mode='constant')
            #H[r*Lg:(r+1)*Lg,:L] = toeplitz(col, row)
            H[r*Lg:(r+1)*Lg,:L] = pra.convmtx(hx, Lg).T

            # build interferer RIR matrix
            hq = pra.lowPassDirac(i_time[r,:,np.newaxis], i_dmp[r,:,np.newaxis], self.Fs, Lh).sum(axis=0)
            #row = np.pad(hq, ((0,L-len(hq))), mode='constant')
            #col = np.pad(hq[:1], ((0, Lg-1)), mode='constant')
            #H[r*Lg:(r+1)*Lg,L:] = toeplitz(col, row)
            H[r*Lg:(r+1)*Lg,L:] = pra.convmtx(hq, Lg).T

        ones = np.ones((K,1))

        # We first assume the sample are uncorrelated
        #K_nq = np.dot(H, H.T) + R_n
        K_nq = np.dot(H[:,L:], H[:,L:].T) + R_n

        # Compute the TD filters
        K_nq_inv = np.linalg.inv(K_nq)
        C = np.dot(K_nq_inv, As)
        B = np.linalg.inv(np.dot(As.T, C))
        g_temp = np.dot(C, np.dot(B, ones))
        self.filters = g_temp.reshape((self.M,Lg))


'''
We create a new Beamformer class for Rake MVDR in time-domain
'''
class RakeMVDR_TD(pra.Beamformer):

    def computeWeights(self, sources, interferers, R_n, epsilon=1e-2):

        dist_mat = pra.distance(self.R, sources)
        s_time = dist_mat / pra.c
        s_dmp = 1./(4*np.pi*dist_mat)

        dist_mat = pra.distance(self.R, interferers)
        i_time = dist_mat / pra.c
        i_dmp = 1./(4*np.pi*dist_mat)

        offset = np.maximum(s_dmp.max(), i_dmp.max())/(np.pi*self.Fs*epsilon)
        t_min = np.minimum(s_time.min(), i_time.min())
        t_max = np.maximum(s_time.max(), i_time.max())

        s_time -= t_min - offset
        i_time -= t_min - offset
        Lh = int((t_max - t_min + 2*offset)*float(self.Fs))

        if (Lh > self.Lg):
            import warnings
            wng = "Beamforming filters length (%d) are shorter than maximum room impulse response length (%d)." % (self.Lg, Lh)
            warnings.warn(wng, UserWarning)

        # the channel matrix
        Lg = self.Lg
        L = self.Lg + Lh - 1
        H = np.zeros((Lg*self.M, 2*L))

        for r in np.arange(self.M):

            hs = pra.lowPassDirac(s_time[r,:,np.newaxis], s_dmp[r,:,np.newaxis], self.Fs, Lh).sum(axis=0)
            row = np.pad(hs, ((0,L-len(hs))), mode='constant')
            col = np.pad(hs[:1], ((0, Lg-1)), mode='constant')
            H[r*Lg:(r+1)*Lg,0:L] = toeplitz(col, row)

            hi = pra.lowPassDirac(i_time[r,:,np.newaxis], i_dmp[r,:,np.newaxis], self.Fs, Lh).sum(axis=0)
            row = np.pad(hi, ((0,L-len(hi))), mode='constant')
            col = np.pad(hi[:1], ((0, Lg-1)), mode='constant')
            H[r*Lg:(r+1)*Lg,L:2*L] = toeplitz(col, row)

        # the constraint vector
        kappa = int(t_max*self.Fs)
        h = H[:,kappa]

        # We first assume the sample are uncorrelated
        Ryy = np.dot(H, H.T) + R_n

        # Compute the TD filters
        Ryy_inv = np.linalg.inv(Ryy)
        g_temp = np.dot(Ryy_inv, h)
        g = g_temp/np.inner(h, g_temp)
        self.filters = g.reshape((self.M,Lg))

        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(np.arange(L)/float(self.Fs), np.dot(H[:,:L].T, g))


