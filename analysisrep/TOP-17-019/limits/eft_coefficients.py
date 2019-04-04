import numpy as np
# from mpl_toolkits.mplot3d import axes3d
# import matplotlib.pyplot as plt
# from matplotlib import cm
# import matplotlib.patches as mpatches
# import matplotlib.lines as mlines

class EftPredictions:
    def __init__(self,wilson_values,cross_section_values,sigma_sm):
        self.wilson_values = wilson_values                  # values of C_i
        self.cross_section_values = cross_section_values    # result of MG calculations for given set of C_i, e.g. [1.,1.,0.,0.,0.]
        self.sig_sm = sigma_sm

        self.sig_i = []     # sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij
        self.sig_ij = []    # sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij
        self.S = []         # sigma_i and sigma_ij in the same vector

        self.N_wilsons = 5
        self.Sigma1_vec  = None
        self.Sigma2_matr = None
        self.name_index_map = {"O_R":0,"O_L^1":1,"O_L^8":2,"O_B^1":3,"O_B^8":4} #indices of Sigma1_vec and Sigma2_matr run in range 0..N_wilsons

        self.calculate_sigmas()

    def gen_row(self,c):
        '''
        Generate single row of the matrix for sigma_i and sigma_ij
        EFT predictions has the following analytical expression:
            sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij
        '''
        row = [c[0], c[1], c[2], c[3], c[4], c[0] ** 2., c[1] ** 2., c[2] ** 2., c[3] ** 2., c[4] ** 2.,
               2. * c[0] * c[1], 2. * c[0] * c[2], 2. * c[0] * c[3], 2. * c[0] * c[4], 2. * c[1] * c[2],
               2. * c[1] * c[3], 2. * c[1] * c[4], 2. * c[2] * c[3], 2. * c[2] * c[4], 2. * c[3] * c[4]]
        return row

    def calculate_sigmas(self):
        # Fill matrix for sigma_i and sigma_ij
        A = np.array( [self.gen_row(c) for c in self.wilson_values] )
        #print A

        # Subtract SM tttt cross section from EFT
        b = []
        for i in range(0, len(self.cross_section_values)): self.cross_section_values[i] = self.cross_section_values[i] * 1000. - self.sig_sm #convert to the same units
        b = np.array(self.cross_section_values)
        #print b

        # solve linear system of equations for sigma_i, sigma_ij coefficients
        S = np.linalg.solve(A, b)
        #print S

        # solution sigma_i, sigma_ij
        self.S = S
        self.sig_i = [S[0], S[1], S[2], S[3], S[4]];
        self.sig_ij = [S[5], S[10], S[11], S[12], S[13], S[10], S[6], S[14], S[15], S[16], S[11], S[14], S[7], S[17], S[18],
                  S[12], S[15], S[17], S[8], S[19], S[13], S[16], S[18], S[19], S[9]];
        #print self.sig_i
        #print self.sig_ij

        self.Sigma1_vec = np.array(self.sig_i)/2.           # the factor of 2 is introduced for later convenience
        self.Sigma2_matr = np.array([[S[5],  S[10], S[11], S[12], S[13]],
                                     [S[10], S[6],  S[14], S[15], S[16]],
                                     [S[11], S[14], S[7],  S[17], S[18]],
                                     [S[12], S[15], S[17], S[8],  S[19]],
                                     [S[13], S[16], S[18], S[19], S[9]] ])


    def CICJ(self,O_R,O_L1,O_L8,O_B1,O_B8):
        '''
        Calculations of the tttt EFT cross section based on just two operators CI, CJ
        :param O_R,O_L1,O_L8,O_B1,O_B8: Wilson coefficients of the corresponding operators
        :param s:
        :return: sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij, where i,j = 0,1
        '''
        c = [O_R,O_L1,O_L8,O_B1,O_B8]
        row = [ self.sig_sm, c[0]*self.S[0], c[1]*self.S[1], c[2]*self.S[2], c[3]*self.S[3], c[4]*self.S[4], 
                (c[0]**2.)*self.S[5], (c[1]**2.)*self.S[6], (c[2]**2.)*self.S[7], (c[3]**2.)*self.S[8], (c[4]**2.)*self.S[9], 
                2.*c[0]*c[1]*self.S[10], 2.*c[0]*c[2]*self.S[11], 2.*c[0]*c[3]*self.S[12], 2.*c[0]*c[4]*self.S[13], 
                2.*c[1]*c[2]*self.S[14], 2.*c[1]*c[3]*self.S[15], 2.*c[1]*c[4]*self.S[16], 
                2.*c[2]*c[3]*self.S[17], 2.*c[2]*c[4]*self.S[18], 
                2.*c[3]*c[4]*self.S[19]]
        return sum(row)

    def O_R_independent_polynomial_coefficients(self):
            '''
            Return a list of polynomial coefficients for 
            tttt cross section parametrisation as function of C0 
            p[0] * C0**n + p[1] * C0**(n-1) + ... + p[n-1]*C0 + p[n]
            '''
            return [self.S[5],self.S[0],self.sig_sm]

    def O_L1_independent_polynomial_coefficients(self):
            '''
            Return a list of polynomial coefficients for 
            tttt cross section parametrisation as function of C1 
            p[0] * C1**n + p[1] * C1**(n-1) + ... + p[n-1]*C1 + p[n]
            '''
            return [self.S[6],self.S[1],self.sig_sm]

    def O_L8_independent_polynomial_coefficients(self):
            '''
            Return a list of polynomial coefficients for 
            tttt cross section parametrisation as function of C2 
            p[0] * C2**n + p[1] * C2**(n-1) + ... + p[n-1]*C2 + p[n]
            '''
            return [self.S[7],self.S[2],self.sig_sm]

    def O_B1_independent_polynomial_coefficients(self):
            '''
            Return a list of polynomial coefficients for 
            tttt cross section parametrisation as function of C3 
            p[0] * C3**n + p[1] * C3**(n-1) + ... + p[n-1]*C3 + p[n]
            '''
            return [self.S[8],self.S[3],self.sig_sm]

    def O_B8_independent_polynomial_coefficients(self):
            '''
            Return a list of polynomial coefficients for 
            tttt cross section parametrisation as function of C4 
            p[0] * C4**n + p[1] * C4**(n-1) + ... + p[n-1]*C4 + p[n]
            '''
            return [self.S[9],self.S[4],self.sig_sm]

    def O_L1L8_independent_polynomial_coefficients(self):
            '''
            Return a list of polynomial coefficients for 
            tttt cross section parametrisation as function of C4 
            p[0] * C4**n + p[1] * C4**(n-1) + ... + p[n-1]*C4 + p[n]
            '''
            return [self.S[6]+self.S[7],self.S[1]+self.S[2],self.sig_sm]      

    def gen_eft_xs(self,wilson_values):
        '''
        Sum individual contributions from different operators

        :param c: vector of Wilson coefficient values, C_i
        :param s: vector of sigma_i, sigma_ij values
        :return:  sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij
        '''
        row = [self.sig_sm, wilson_values[0] * self.S[0], wilson_values[1] * self.S[1], wilson_values[2] * self.S[2], wilson_values[3] * self.S[3], wilson_values[4] * self.S[4],
               (wilson_values[0] ** 2.) * self.S[5], (wilson_values[1] ** 2.) * self.S[6], (wilson_values[2] ** 2.) * self.S[7], (wilson_values[3] ** 2.) * self.S[8], (wilson_values[4] ** 2.) * self.S[9],
               2. * wilson_values[0] * wilson_values[1] * self.S[10], 2. * wilson_values[0] * wilson_values[2] * self.S[11], 2. * wilson_values[0] * wilson_values[3] * self.S[12], 2. * wilson_values[0] * wilson_values[4] * self.S[13],
               2. * wilson_values[1] * wilson_values[2] * self.S[14], 2. * wilson_values[1] * wilson_values[3] * self.S[15], 2. * wilson_values[1] * wilson_values[4] * self.S[16],
               2. * wilson_values[2] * wilson_values[3] * self.S[17], 2. * wilson_values[2] * wilson_values[4] * self.S[18],
               2. * wilson_values[3] * wilson_values[4] * self.S[19]]
        return sum(row)

    def vgen_eft_xs(self,C0,C1,C2,C3,C4):
        '''
        Sum individual contributions from different operators.
        This is the same as gen_eft_xs(), but numpy friendly

        :param c: vector of Wilson coefficient values, C_i
        :param s: vector of sigma_i, sigma_ij values
        :return:  sigma_tttt = sigma_tttt_SM + sum_i C_i*sigma_i + sum_ij C_i*C_j*sigma_ij
        '''
        row = [self.sig_sm, C0 * self.S[0], C1 * self.S[1], C2 * self.S[2], C3 * self.S[3], C4 * self.S[4],
               (C0 ** 2.) * self.S[5], (C1 ** 2.) * self.S[6], (C2 ** 2.) * self.S[7], (C3 ** 2.) * self.S[8], (C4 ** 2.) * self.S[9],
               2. * C0 * C1 * self.S[10], 2. * C0 * C2 * self.S[11], 2. * C0 * C3 * self.S[12], 2. * C0 * C4 * self.S[13],
               2. * C1 * C2 * self.S[14], 2. * C1 * C3 * self.S[15], 2. * C1 * C4 * self.S[16],
               2. * C2 * C3 * self.S[17], 2. * C2 * C4 * self.S[18],
               2. * C3 * C4 * self.S[19]]
        return sum(row)
