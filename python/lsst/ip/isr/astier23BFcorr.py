import numpy as np
import pyfftw


# code adaped from :
# https://stackoverflow.com/questions/14786920/convolution-of-two-three-dimensional-arrays-with-padding-on-one-side-too-slow
# code posted by Henry Gomersal 
class CustomFFTConvolution(object):
    """
    A class that performs image convolutions in Fourier space, using pyfftw.
    The constructor takes images as arguments, and creates the 
    plans in fftw3 parlance. The convolutions are done by the __call_ routine.
    This is faster than scipy.signal.fftconvole, and it saves some transforms
    by allowing to convolve the same image with several kernels.
    pyfftw does not accomodate float32 images, so everything 
    should be double precision.
    """
    def __init__(self, A, B, threads=1):
        # minimum size of the convolution
        shape = (np.array(A.shape) + np.array(B.shape))-1
        # immediate larger "fast size". Can be a huge gain.
        shape = np.array([pyfftw.next_fast_len(s) for s in shape])
        # fftw cooks up plans:
        self.fft_plan_im = pyfftw.builders.rfftn(
                    A, s=shape, threads=threads)
        self.fft_plan_kern = pyfftw.builders.rfftn(
                    B, s=shape, threads=threads)
        self.ifft_plan = pyfftw.builders.irfftn(
                    self.fft_plan_im.get_output_array(), s=shape,
                    threads=threads)

    def __call__(self, im, kernels):
        """
        Carries out the convolution and trims the result to the size of im.
        if kernels is a list, then the routine returns 
        the list of corresponding convolutions.
        """
        # accomodate both a list of kernels and a single kernel
        l = [kernels] if type(kernels) != list else kernels
        convs = []
        for kern in l:
            # transform the image and the kernel
            tim = self.fft_plan_im(im)
            tkern = self.fft_plan_kern(kern)

            conv = self.ifft_plan(tim*tkern)
            # now trim the result
            # follow the 'same' policy of scipy.signal.fftconvolve
            oy = kern.shape[0]//2
            ox = kern.shape[1]//2
            convs.append(conv[oy:oy+im.shape[0], ox:ox+im.shape[1]].copy())
        return convs[0] if type(kernels) != list else convs


class BFCorr:
    """
    Evaluates the correction of CCD images affected by the
    brighter-fatter effect, along what is described in
    https://arxiv.org/abs/2301.03274. Requires as input the result of
    an electrostatic fit to flat covariance data (or any other
    determination of pixel boundary shifts under the influcence of a
    single electron)
    """
    # by discussing with Pierre Astier, last version of the files are
    # located here at s3df: /sdf/home/a/astier/place/run7/E2016/R??_S??/avalues.npy
    def __init__(self, filename):
        """
        Filename  refers to an input tuple that contains the 
        boundary shifts for one electron. This file is produced by an
        electrostatic fit to data extracted from flat-field statistics,
        implemented in https://gitlab.in2p3.fr/astier/bfptc/tools/fit_cov.py
        """
        n = np.load(filename).view(np.recarray)
        # n = croaks.NTuple.fromtxt(filename).view(np.recarray)
        self.range = np.max(n.i)
        r = self.range
        # the input file contains one quadrant and we need 4 for convolutions.
        self.kN = np.zeros((2*r+1,2*r+1))
        self.kE = np.zeros_like(self.kN)
        # self.kW = np.zeros_like(self.kN)
        # self.kS = np.zeros_like(self.kN)

        # fill the 4 quadrants
        # i refers to serial direction, j to parallel
        self.kN[r + n.i, r + n.j] = n.aN
        self.kN[r - n.i, r + n.j] = n.aN
        self.kN[r + n.i, r - n.j] = n.aS
        self.kN[r - n.i, r - n.j] = n.aS
        self.kE[r + n.i, r + n.j] = n.aE
        self.kE[r + n.i, r - n.j] = n.aE
        self.kE[r - n.i, r + n.j] = n.aW
        self.kE[r - n.i, r - n.j] = n.aW
        # tweak the edges so that the sum rule applies.
        self.kN[:, 0] = -self.kN[:,-1]
        self.kE[0, :] = -self.kE[-1,:]
        print("INFO: BF kernel sum rules : kN %f, kE %f"%(self.kN.sum(), self.kE.sum()))
                
        # We use the normalization of Guyonnet et al (2015)
        # (compatible with the way the input file is produced).
        # 1/2 is due to the fact that the charge distribution at the end
        # is twice of the average, and the second 1/2 is due to
        # charge interpolation.
        self.kN *= 0.25
        self.kE *= 0.25
        # indeed, i and j in the tuple refer to serial and parallel directions
        # in most of the python codes, the imeage reads im[j,i], so :
        self.kN = self.kN.T
        self.kE = self.kE.T

    
    def DeltaImageFFT(self, im):
        """
        Computes the correction and returns the "delta_image",
        to be subtracted from "im" in order to *undo* the BF effect.
        im should be expressed in *electrons*, and hence can contain
        several segments corresponding to different video channels.
        The returned image is also expressed in electrons.
        im is unchanged.
        """
        im = im.astype('float64') # mandatory for fftw
        conv_obj = CustomFFTConvolution(im, self.kN)
        convs = conv_obj(im, [self.kN, self.kE])
        # convs contains the boundary shifts (in pixels size units)
        # for [horizontal, vertical] boundaries.
        # we now compute the charge to move around 
        delta = np.zeros_like(im)
        boundary_charge = np.zeros_like(im)

        # horizontal boundaries (// direction)
        # we could use a more elaborate interpolator for estimating the
        # charge on the boundary
        boundary_charge[:-1,:] = im[1:,:]+im[:-1,:]
        # boundary_charge[1:-2,:] = (9./8.)*(I[2:-1,:]+I[1:-2,:] -
        # (1./8.)*(I[0:-3,:]+I[3:,:])
        
        # the charge to move around is the
        # product of the boundary shift (in pixel size unit) times the
        # charge on the boundary (in charge per pixel unit)
        dq = boundary_charge*convs[0]
        delta += dq
        # what is gained by a pixel is lost by its neighbor (the righ one!)
        delta[1:,:] -= dq[:-1,:]

        # vertical boundaries
        boundary_charge = np.zeros_like(im) #  reset to zero
        # same comment as above
        boundary_charge[:,:-1] = im[:,1:]+im[:,:-1]
        dq = boundary_charge*convs[1]
        delta += dq
        # what is gained by a pixel is lost by its neighbor
        delta[:,1:] -=  dq[:,:-1]
        # one might check that delta.sum() ~ 0 (charge conservation)
        return delta
