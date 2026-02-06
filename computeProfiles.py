from headers import *

##################################################################################
##################################################################################


# def computeProfiles(ts, filterType, est='tsz_uniformweight', iBootstrap=None, iVShuffle=None, tTh='',  mVir=None, z=[0., 100.], mask=None):
#     """Returns the estimated profile and its uncertainty for each aperture.
#     est: string to select the estimator
#     iBootstrap: index for bootstrap resampling
#     iVShuffle: index for shuffling velocities
#     tTh: to replace measured temperatures by a theory expectation
#     ts: option to specify another thumbstack object
#     """

#     # tStart = time()

#     print(("- Compute stacked profile: "+filterType+", "+est+", "+tTh))

#     # compute stacked profile from another thumbstack object
#     if mVir is None:
#         mVir = [ts.mMin, ts.mMax]

#     # select objects that overlap, and reject point sources
#     if mask is None:
#         mask = ts.catalogMask(overlap=True, psMask=True,
#                               filterType=filterType, mVir=mVir, z=z)

# #      tMean = ts.meanT[filterType].copy()

#     # temperatures [muK * sr]
#     if tTh == '':
#         t = ts.filtMap[filterType].copy()  # [muK * sr]
#     elif tTh == 'tsz':
#         # expected tSZ signal
#         # AP profile shape, between 0 and 1
#         sigma_cluster = 3.  # 1.5  # arcmin
#         shape = ts.ftheoryGaussianProfile(
#             sigma_cluster)  # between 0 and 1 [dimless]
#         # multiply by integrated y to get y profile [sr]
#         t = np.column_stack(
#             [ts.Catalog.integratedY[:] * shape[iAp] for iAp in range(ts.nRAp)])
#         # convert from y profile to dT profile if needed
#         if ts.cmbUnitLatex == r'$\mu$K':
#             nu = ts.cmbNu   # Hz
#             Tcmb = 2.726   # K
#             h = 6.63e-34   # SI
#             kB = 1.38e-23  # SI

#             def f(nu):
#                 """frequency dependence for tSZ temperature
#                 """
#                 x = h*nu/(kB*Tcmb)
#                 return x*(np.exp(x)+1.)/(np.exp(x)-1.) - 4.
#             # t *= 2. * f(nu) * Tcmb * 1.e6  # [muK * sr]
#             t *= f(nu) * Tcmb * 1.e6  # [muK * sr]
#     t = t[mask, :]
# #     tMean = tMean[mask,:]
#     # -v/c [dimless]
#     v = -ts.Catalog.vR[mask] / 3.e5
#     v -= np.mean(v)

#     # true filter variance for each object and aperture,
#     # valid whether or not a hit count map is available
#     s2Full = ts.filtVarTrue[filterType][mask, :]
#     # Variance from hit count (if available)
#     s2Hit = ts.filtHitNoiseStdDev[filterType][mask, :]**2
#     # print "Shape of s2Hit = ", s2Hit.shape
#     # halo masses
#     m = ts.Catalog.Mvir[mask]

#     if iBootstrap is not None:
#         # make sure each resample is independent,
#         # and make the resampling reproducible
#         np.random.seed(iBootstrap)
#         # list of overlapping objects
#         nObj = np.sum(mask)
#         # print "sample "iBootstrap, ";", nObj, "objects overlap with", ts.name
#         I = np.arange(nObj)
#         # choose with replacement from this list
#         J = np.random.choice(I, size=nObj, replace=True)
#         #
#         t = t[J, :]
#         # tMean = tMean[J,:]
#         v = v[J]
#         s2Hit = s2Hit[J, :]
#         s2Full = s2Full[J, :]
#         m = m[J]

#     if iVShuffle is not None:
#         # make sure each shuffling is independent,
#         # and make the shuffling reproducible
#         np.random.seed(iVShuffle)
#         # list of overlapping objects
#         nObj = np.sum(mask)
#         I = np.arange(nObj)
#         # shuffle the velocities
#         J = np.random.permutation(I)
#         #
#         v = v[J]

#     # tSZ: uniform weighting
#     if est == 'tsz_uniformweight':
#         weights = np.ones_like(s2Hit)
#         norm = 1./np.sum(weights, axis=0)
#     # tSZ: detector-noise weighted (hit count)
#     elif est == 'tsz_hitweight':
#         weights = 1./s2Hit
#         norm = 1./np.sum(weights, axis=0)
#     # tSZ: full noise weighted (detector noise + CMB)
#     elif est == 'tsz_varweight':
#         weights = 1./s2Full
#         norm = 1./np.sum(weights, axis=0)

#     # kSZ: velocity-weighted uniform weighting
    


#     # tStop = time()
#     # print "stacked profile took", tStop-tStart, "sec"

#     # return the stacked profiles
#     # stack = norm * np.sum(t * weights, axis=0)
#     # sStack = norm * np.sqrt(np.sum(s2Full * weights**2, axis=0))
#     return t

def computeProfiles(ts, filterType, est, iBootstrap=None, iVShuffle=None, tTh='', mVir=None, z=[0., 100.], mask=None):
    """Returns the estimated profile and its uncertainty for each aperture.
    est: string to select the estimator
    iBootstrap: index for bootstrap resampling
    iVShuffle: index for shuffling velocities
    tTh: to replace measured temperatures by a theory expectation
    ts: option to specify another thumbstack object
    """

    #tStart = time()

    #print(("- Compute stacked profile: "+filterType+", "+est+", "+tTh))

    # compute stacked profile from another thumbstack object
    # if ts is None:
        # ts = self
    if mVir is None:
        mVir = [ts.mMin, ts.mMax]

    # select objects that overlap, and reject point sources
    if mask is None:
        mask = ts.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=mVir, z=z)

    #      tMean = ts.meanT[filterType].copy()

    # temperatures [muK * sr]
    if tTh=='':
        t = ts.filtMap[filterType].copy() # [muK * sr]
    elif tTh=='tsz':
        # expected tSZ signal
        # AP profile shape, between 0 and 1
        sigma_cluster = 3.   #1.5  # arcmin
        shape = ts.ftheoryGaussianProfile(sigma_cluster) # between 0 and 1 [dimless]
        # multiply by integrated y to get y profile [sr]
        t = np.column_stack([ts.Catalog.integratedY[:] * shape[iAp] for iAp in range(ts.nRAp)])
        # convert from y profile to dT profile if needed
        if ts.cmbUnitLatex==r'$\mu$K':
            nu = ts.cmbNu   # Hz
            Tcmb = 2.726   # K
            h = 6.63e-34   # SI
            kB = 1.38e-23  # SI
            def f(nu):
                """frequency dependence for tSZ temperature
                """
                x = h*nu/(kB*Tcmb)
                return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
            #t *= 2. * f(nu) * Tcmb * 1.e6  # [muK * sr]
            t *= f(nu) * Tcmb * 1.e6  # [muK * sr]
    elif tTh=='ksz':
        # expected kSZ signal
        # AP profile shape, between 0 and 1
        sigma_cluster = 1.5  # arcmin
        shape = ts.ftheoryGaussianProfile(sigma_cluster) # between 0 and 1 [dimless]
        # multiply by integrated kSZ to get kSZ profile [muK * sr]
        t = np.column_stack([ts.Catalog.integratedKSZ[:] * shape[iAp] for iAp in range(ts.nRAp)])   # [muK * sr]
        if ts.cmbUnitLatex=='':
            t /= 2.726e6   # convert from [muK*sr] to [sr]
    t = t[mask, :]
    #     tMean = tMean[mask,:]
    # -v/c [dimless]
    v = -ts.Catalog.vR[mask] / 3.e5
    #v -= np.mean(v) # tuks!!!!!!!!!!!!!!!!!!! TESTING
    # I don't think we should do this!!!!!!!!!!!!! tuks B.H.

    #      # expected sigma_{v_{true}}, for the normalization
    #      #print "computing v1d norm"
    #      #tStartV = time()
    #      z = ts.Catalog.Z[mask]
    #      #f = lambda zGal: ts.U.v1dRms(0., zGal, W3d_sth)**2
    #      #sVTrue = np.sqrt(np.mean(np.array(map(f, z))))
    #      sVTrue = ts.U.v1dRms(0., np.mean(z), W3d_sth) / 3.e5  # (v^true_rms/c) [dimless]
    #      #print "sigma_v_true =", sVTrue
    #      #print "at z=0.57, expect", np.sqrt(f(0.57))
    #      #tStopV = time()
    #      #print "v1d norm took", tStopV - tStartV, "sec"

    #true filter variance for each object and aperture,
    # valid whether or not a hit count map is available
    s2Full = ts.filtVarTrue[filterType][mask, :]
    # Variance from hit count (if available)
    s2Hit = ts.filtHitNoiseStdDev[filterType][mask, :]**2
    #print "Shape of s2Hit = ", s2Hit.shape
    # halo masses
    m = ts.Catalog.Mvir[mask]

    if iBootstrap is not None:
        # make sure each resample is independent,
        # and make the resampling reproducible
        np.random.seed(iBootstrap)
        # list of overlapping objects
        nObj = np.sum(mask)
        #print "sample "iBootstrap, ";", nObj, "objects overlap with", ts.name
        I = np.arange(nObj)
        # choose with replacement from this list
        J = np.random.choice(I, size=nObj, replace=True)
        #
        t = t[J,:]
        #tMean = tMean[J,:]
        v = v[J]
        s2Hit = s2Hit[J,:]
        s2Full = s2Full[J,:]
        m = m[J]

    if iVShuffle is not None:
        # make sure each shuffling is independent,
        # and make the shuffling reproducible
        np.random.seed(iVShuffle)
        # list of overlapping objects
        nObj = np.sum(mask)
        I = np.arange(nObj)
        # shuffle the velocities
        J = np.random.permutation(I)
        #
        v = v[J]

    # tSZ: uniform weighting
    if est=='tsz_uniformweight' or est=='tsz_anisotropic_uniformweight':
        weights = np.ones_like(s2Hit)
        norm = 1./np.sum(weights, axis=0)
    # tSZ: detector-noise weighted (hit count)
    elif est=='tsz_hitweight' or est=='tsz_anisotropic_hitweight':
        weights = 1./s2Hit
        norm = 1./np.sum(weights, axis=0)
    # tSZ: full noise weighted (detector noise + CMB)
    elif est=='tsz_varweight' or est=='tsz_anisotropic_varweight':
        weights = 1./s2Full
        norm = 1./np.sum(weights, axis=0)

    # tau: uniform weighting
    elif est=='tau_uniformweight' or est=='tau_anisotropic_uniformweight':
        weights = np.sign(v[:, np.newaxis]) * np.ones_like(s2Hit) # S[T_l]
        norm = 1./(np.mean(np.abs(v)) * np.sum(np.abs(weights), axis=0)) # 1/(<|T_l|> N)

    # kSZ: uniform weighting
    elif est=='ksz_uniformweight' or est=='ksz_anisotropic_uniformweight':
        # remove mean temperature
        #t -= np.mean(t, axis=0)
        #         t -= tMean
        weights = v[:,np.newaxis] * np.ones_like(s2Hit)
        #norm = sVTrue / np.sum(v[:,np.newaxis]*weights, axis=0)
        norm = np.std(v) / ts.Catalog.rV / np.sum(v[:,np.newaxis]*weights, axis=0)
    # kSZ: detector-noise weighted (hit count)
    elif est=='ksz_hitweight' or est=='ksz_anisotropic_hitweight':
        # remove mean temperature
        #t -= np.mean(t, axis=0)
        #         t -= tMean
        weights = v[:,np.newaxis] / s2Hit
        #norm = sVTrue / np.sum(v[:,np.newaxis]*weights, axis=0)
        norm = np.std(v) / ts.Catalog.rV / np.sum(v[:,np.newaxis]*weights, axis=0)
    # kSZ: full noise weighted (detector noise + CMB)
    elif est=='ksz_varweight' or est=='ksz_anisotropic_varweight':
        # remove mean temperature
        #t -= np.mean(t, axis=0)
        #         t -= tMean
        weights = v[:,np.newaxis] / s2Full
        #norm = sVTrue / np.sum(v[:,np.newaxis]*weights, axis=0)
        norm = np.std(v) / ts.Catalog.rV / np.sum(v[:,np.newaxis]*weights, axis=0)
    # kSZ: full noise weighted (detector noise + CMB)
    elif est=='ksz_massvarweight' or est=='ksz_anisotropic_massvarweight':
        # remove mean temperature
        #t -= np.mean(t, axis=0)
        #         t -= tMean
        weights = m[:,np.newaxis] * v[:,np.newaxis] / s2Full
        #norm = np.mean(m) * sVTrue / np.sum(m[:,np.newaxis]**2 * v[:,np.newaxis]**2 / s2Full, axis=0)
        norm = np.mean(m) * np.std(v) / ts.Catalog.rV / np.sum(m[:,np.newaxis]**2 * v[:,np.newaxis]**2 / s2Full, axis=0)

    #tStop = time()
    #print "stacked profile took", tStop-tStart, "sec"

    # return all the profiles
    # # stack = norm * np.sum(t * weights, axis=0)
    # sStack = norm * np.sqrt(np.sum(s2Full * weights**2, axis=0))
    return t * weights

    '''
    # or, if requested, compute and return the stacked cutout map
    else:
        # define chunks
        nChunk = ts.nProc
        chunkSize = int(ts.Catalog.nObj / nChunk)
        # list of indices for each of the nChunk chunks
        chunkIndices = [list(range(iChunk*chunkSize, (iChunk+1)*chunkSize)) for iChunk in range(nChunk)]
        # make sure not to miss the last few objects:
        # add them to the last chunk
        chunkIndices[-1] = list(range((nChunk-1)*chunkSize, ts.Catalog.nObj))

        # select weights for a typical aperture size (not the smallest, not the largest)
        #iRAp0 = ts.nRAp / 2
        iRAp0 = int(ts.nRAp / 4)
        norm = norm[iRAp0]
        # need to link object number with weight,
        # despite the mask
        weightsLong = np.zeros(ts.Catalog.nObj)
        weightsLong[mask] = weights[:,iRAp0]

        def stackChunk(iChunk):
            # object indices to be processed
            chunk = chunkIndices[iChunk]

            # start with a null map for stacking
            resMap = ts.cutoutGeometry()

            # radian positions of each pixel
            ipos = resMap.posmap()
            X = ipos[0]
            Y = ipos[1]
            x, y = X[:, 0], Y[0, :]
            XY = np.array([X.flatten(), Y.flatten()])

            # TESTING!!!!!!!!!!!!!! I am feeling a bit lazy tbqh....
            # size of canvas in radians (this is for just x direction, but cutout is symmetric)
            size = ipos[0,:,:].max() - ipos[0,:,:].min()

            # size of pixel in radians
            dx = float(size) / (resMap.shape[0]-1)
            dy = float(size) / (resMap.shape[1]-1)

            # centers of the pixels (index+0.5 times size)
            x_bins = dx * (np.arange(resMap.shape[0]+1) - 0.5)
            y_bins = dy * (np.arange(resMap.shape[1]+1) - 0.5)
            x_grid, y_grid = np.meshgrid(x_bins, y_bins, indexing='ij')
            #XY = np.array([X.flatten(), Y.flatten()]) # same? I actually don't think this is right

            x_grid = x_grid*180./np.pi*60.
            y_grid = y_grid*180./np.pi*60.
            cell_size = x_grid[1, 0]-x_grid[0, 0]
            x_grid += cell_size/2.
            y_grid += cell_size/2.
            x_grid = x_grid[:-1, :-1]
            y_grid = y_grid[:-1, :-1]
            x_grid -= (x_grid.max()-x_grid.min())/2.
            y_grid -= (y_grid.max()-y_grid.min())/2.

            r = np.sqrt(x_grid**2+y_grid**2)
            #mu = np.abs(y_grid/r)
            th = np.arctan2(y_grid, x_grid)
            th[th < 0.] += 2.*np.pi

            """
            from scipy.special import legendre
            order = 2
            Ln = legendre(order)
            """

            r_bins = np.linspace(1., np.floor(x_grid.max()), 11)#6)
            r_binc = (r_bins[:-1] + r_bins[1:]) / (2.0)
            """
            mu_bins = np.linspace(0., 1., 10)
            mu_binc = (mu_bins[:-1] + mu_bins[1:]) / (2.0)
            hist_norm, _, _ = np.histogram2d(r.flatten(), mu.flatten(), bins=[r_bins, mu_bins])
            """
            hist_norm, _ = np.histogram(r.flatten(), bins=r_bins)

            m_ell0 = np.zeros((len(chunk), len(r_binc)))
            m_ell2 = np.zeros((len(chunk), len(r_binc)))
            count = 0
            want_random = False # TESTING
            if want_random:
                seed_def = 6000 # def
                seed = 3000
                np.random.seed(seed) # randomized
            for iObj in chunk:
                if iObj%10000==0:
                    print("- analyze object", iObj)
                if ts.overlapFlag[iObj]:
                    # Object coordinates
                    ra = ts.Catalog.RA[iObj]   # in deg
                    dec = ts.Catalog.DEC[iObj] # in deg
                    z = ts.Catalog.Z[iObj] # not used
                    # extract postage stamp around it
                    opos, stampMap, stampMask, stampHit = ts.extractStamp(ra, dec, test=False)
                    if "anisotropic" in est:
                        # randomized
                        if want_random:
                            alpha = np.random.rand()*2.*np.pi
                            ca = np.cos(alpha)
                            sa = np.sin(alpha)
                        else:
                            ca = ts.Catalog.vX[iObj] # cos(alpha)
                            sa = ts.Catalog.vY[iObj] # sin(alpha)
                        fun2D = RectBivariateSpline(x, y, stampMap, kx=1, ky=1)
                        R = np.array([[ca, sa], [-sa, ca]]) # tested that this is the right
                        #R = np.array([[ca, -sa], [sa, ca]]) # TESTING!!!! I think mirror reflected
                        X_rot, Y_rot = np.dot(R, XY)
                        stampMap = fun2D(X_rot, Y_rot, grid=False).reshape(resMap.shape)
                        del X_rot, Y_rot
                    resMap += (stampMap * weightsLong[iObj])

                    if "anisotropic" in est:
                        m0, _ = np.histogram(r.flatten(), bins=r_bins, weights=(stampMap * weightsLong[iObj]).flatten())
                        #m0, _ = np.histogram(r.flatten(), bins=r_bins, weights=(stampMap).flatten()) # now
                        #m0, _ = np.histogram(r.flatten(), bins=r_bins, weights=(stampMap*np.sign(np.random.rand()-0.5)).flatten()) # randw
                        m2, _ = np.histogram(r.flatten(), bins=r_bins, weights=(stampMap * weightsLong[iObj]).flatten()*np.cos(2.*th.flatten()))
                        #m2, _ = np.histogram(r.flatten(), bins=r_bins, weights=(stampMap).flatten()*np.cos(2.*th.flatten())) # now
                        #m2, _ = np.histogram(r.flatten(), bins=r_bins, weights=(stampMap*np.sign(np.random.rand()-0.5)).flatten()*np.cos(2.*th.flatten())) # randw
                        #m4, _ = np.histogram(r.flatten(), bins=r_bins, weights=(stampMap * weightsLong[iObj]).flatten()*np.cos(4.*th.flatten()))
                        m0 /= hist_norm
                        m2 /= hist_norm/2.
                        #m4 /= hist_norm/2.
                        m_ell0[count] = m0
                        m_ell2[count] = m2

                count += 1
            if "anisotropic" in est:
                # save the m_ell
                # randomized
                if want_random:
                    random_str = "_random"
                    if seed != seed_def: random_str += f"{seed}" # type: ignore
                else:
                    random_str = ""
                #smooth_str = "_r006"
                smooth_str = "_r003"
                if ts.cmbNu == 150.e9:
                    version_str = "_dr5_f150"
                elif ts.cmbNu == 90.e9:
                    version_str = "_dr5"
                else:
                    version_str = "_dr5"
                if "varweight" in est:
                    version_str += "_varweight"
                #version_str += "_now"
                #version_str += "_randw"
                mell_dir = f"/pscratch/sd/b/boryanah/ACTxDESI/anisotropic{random_str}_wide{smooth_str}{version_str}/"
                os.makedirs(mell_dir, exist_ok=True)
                np.savez(f"{mell_dir}/mell2_{iChunk:d}.npz", m_ell0=m_ell0, m_ell2=m_ell2, r_bins=r_bins, gal_inds=chunk)
            return resMap

        # dispatch each chunk of objects to a different processor
        with sharedmem.MapReduce(np=ts.nProc) as pool:
            resMap = np.array(pool.map(stackChunk, list(range(nChunk))))

        # sum all the chunks
        resMap = np.sum(resMap, axis=0)
        # normalize by the proper sum of weights
        resMap *= norm

        return resMap
    '''