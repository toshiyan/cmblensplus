
# Define quad estimator names
class quad_fname:

    def __init__(self,qtype,qest,root,qtag,otag,ids):

        # qtype is the type of mode coupling, such as lens, rot, etc
        qalm = root + qtype + '/alm/'
        qrdn = root + qtype + '/rdn0/'
        qmlm = root + qtype + '/mean/'
        qaps = root + qtype + '/aps/'

        # normalization and tau transfer function
        self.al   = qaps+'Al_'+qest+'_'+qtag+'.dat'
        self.wl   = qaps+'Wl_'+qest+'_'+qtag+'.dat'

        # N0/N1 bias
        self.n0bl = qaps+'n0_'+qest+'_'+qtag+'.dat'
        self.n1bs = qaps+'n1_'+qest+'_'+qtag+'.dat'

        # mean field
        self.ml   = [qmlm+'cl_'+qest+'_'+qtag+'_'+x+'.dat' for x in ids]
        self.mfb  = [qmlm+'mlm_'+qest+'_'+qtag+'_'+x+'.pkl' for x in ids]

        # reconstructed spectra
        self.mcls = qaps+'cl_'+qest+'_'+qtag+'.dat'
        self.mcbs = qaps+'cl_'+qest+'_'+qtag+otag+'.dat'
        self.ocls = qaps+'cl_'+ids[0]+'_'+qest+'_'+qtag+'.dat'
        self.ocbs = qaps+'cl_'+ids[0]+'_'+qest+'_'+qtag+otag+'.dat'
        self.cl   = [qaps+'rlz/cl_'+qest+'_'+qtag+'_'+x+'.dat' for x in ids]

        # reconstructed alm and RDN0
        self.alm  = [qalm+'alm_'+qest+'_'+qtag+'_'+x+'.pkl' for x in ids]
        self.rdn0 = [qrdn+'rdn0_'+qest+'_'+qtag+'_'+x+'.dat' for x in ids]



