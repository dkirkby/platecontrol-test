import numpy as np
import posconstants as pc


# %% rotation matrices
def Rx(angle):  # all in radians
    Rx = np.array([
                   [1.0,           0.0,            0.0],
                   [0.0,           np.cos(angle),  -np.sin(angle)],
                   [0.0,           np.sin(angle),  np.cos(angle)]
                  ])
    return Rx


def Ry(angle):  # all in radians
    Ry = np.array([
                   [np.cos(angle),  0.0,            np.sin(angle)],
                   [0.0,            1.0,            0.0],
                   [-np.sin(angle), 0.0,            np.cos(angle)]
                  ])
    return Ry


def Rz(angle):  # all in radians
    Rz = np.array([
                   [np.cos(angle), -np.sin(angle), 0.0],
                   [np.sin(angle), np.cos(angle),  0.0],
                   [0.0,           0.0,            1.0]
                  ])
    return Rz


def Rxyz(alpha, beta, gamma):  # yaw-pitch-roll system, all in radians
    return Rz(gamma) @ Ry(beta) @ Rx(alpha)  # @ is matrix multiplication


def typecast(x):
    if type(x) in [list, tuple]:  # 2 or 3 coordinates, reshape into column vec
        return np.array(x).reshape(len(x), 1)
    elif type(x) is np.ndarray:
        if x.ndim == 1:  # a 1D array of 2 or 3 coordinates
            return x.reshape(x.size, 1)  # reshape into column vector
        elif x.ndim == 2:  # 2D array, assume each column is a column vec
            return x
        else:
            raise ValueError(f'Wrong numpy array dimension: {x.ndim}')
    else:
        raise TypeError(f'Wrong data type: {type(x)}')


# %% class definition
class PetalTransforms:
    """
    This class provides transformations between the petal local nominal
    coordinate system (CAD model) and the focal plane CS5.
    The transformation parameters are found in focal palte alignment as
    6 DOF rigid-body transformation parameters, defined in trans_keys below,
    and reported per petal according to DESI-2850.

    trans_keys = ['petal_offset_x', 'petal_offset_y', 'petal_offset_z',
                  'petal_rot_x', 'petal_rot_y', 'petal_rot_z']
    which are the input variables, Tx, Ty, Tz in mm, alpha, beta, gamma in rad

    For forward transformation,
    the transformation defined by these 6 parameters can act on
    (1) nominal values in the petal local CS as defined in the CAD model, or
    (2) metrology data in the petal local CS (best-fit to the CAD, "ZBF")
    and result in coordinates in CS5, assuming the alignment achieved on
    mountain is the same as done in the metrology lab.

    In use, you make one instance of PetalTransforms per petal.
    A given instance can be initialised with either 6 dof or just 3 dof as
    viewed by the FVC, assuming the other three are zero.

    Note the rotation angles reported from metrology and used for
    focal plate alingment are NOT the classic Euler 3-2-3 or Z-Y-Z angles.
    They are the rotation angles w.r.t. x, y, and z-axes.
    The Euler method provided in this class contains the
    canonical order of operations to do this transform in 3-2-3.

    The coordinate systems are:

        ptlXYZ:     position of a point in the petal's local coordinate system
                    according to the CAD design (overlaps with location 3)

        obsXYZ:     (x, y, z) global cartesian coordinates on the focal plate
                    cenetered on optical axis looking at fiber tips (don't use)

        QST:        (q,s, t) global coordinates on the focal plate, per
                    DESI-0742, 5288
                    q is the angle about the optical axis
                    s is the path distance from optical axis, along the curved
                        focal surface,
                    t is the local normal vector out of the focal surface

        flatXY:     reserved for petal local CS, with focal plate slightly
                    stretched out and flattened (used for anti-collision)

                    for global flatXY coordinates, use obsXY methods with
                    curved = False

    With QS coordinates, PosTransforms converts into any other
    positioner-specific CS. The wrappers here convert a list of positioner
    coordinates, but that is not acturally supported in PosTransforms.
    We may want to get rid them and unify our transformation methods.

    All methods take list or np.ndarray input as handled in typecast() above.
    Set optional arg cast=True when feeding a non-column-vector-based array.

    The option "curved" sets whether we are looking at a flat focal plane
    (such as a laboratory test stand), or a true petal, with its asphere.
    When curved == False, PosTransforms will treat S as exactly equal to the
    radius R.
    """

    dtb_device_ids = [543, 544, 545]  # three datum tooling balls on a petal

    def __init__(self, Tx=0, Ty=0, Tz=0, alpha=0, beta=0, gamma=0,
                 curved=True):  # input in mm and radians
        # self.postrans = PosTransforms(curved=True)
        # tanslation matrix from petal nominal CS to CS5, column vector
        self.T = np.array([Tx, Ty, Tz]).reshape(3, 1)
        # orthogonal rotation matrix
        self.R = Rxyz(alpha, beta, gamma)
        self.curved = curved
        self.petal_alignment = {'Tx': Tx, 'Ty': Ty, 'Tz': Tz,
                                'alpha': alpha, 'beta': beta, 'gamma': gamma}

    # %% QST transforms in 3D
    @staticmethod
    def QST_to_obsXYZ(QST, cast=False):
        """Transforms QST coordinates into obsXYZ system.
        INPUT:  3 x N array, each column vector is QST
        OUTPUT: 3 x N array, each column vector is obsXYZ
        """
        if cast:
            QST = typecast(QST)
        Q_rad, S, T = np.radians(QST[0, :]), QST[1, :], QST[2, :]
        R = pc.S2R_lookup(S) + T * np.sin(np.radians(pc.S2N_lookup(S)))
        Z = pc.S2Z_lookup(S) + T * np.cos(np.radians(pc.S2N_lookup(S)))
        X, Y = R * np.cos(Q_rad), R * np.sin(Q_rad)
        return np.vstack([X, Y, Z])  # 3 x N array

    @staticmethod
    def obsXYZ_to_QST(obsXYZ, cast=False, max_iter=10, tol=1e-6):
        """Transform obsXYZ coordinates into QST system, tol of dtheta in deg
        INPUT:  3 x N array, each column vector is obsXYZ
        OUTPUT: 3 x N array, each column vector is QST
        """
        if cast:
            obsXYZ = typecast(obsXYZ)
        X, Y, Z = obsXYZ[0, :], obsXYZ[1, :], obsXYZ[2, :]
        Q = np.degrees(np.arctan2(Y, X))  # Y over X, azimuthal angle in deg
        R = np.sqrt(np.square(X) + np.square(Y))
        # iteratively find S, use nutation which is the fastest changing coord
        N0 = pc.R2N_lookup(R)  # get nutation angle in the ballpark
        for i in range(max_iter):
            # print(f'iter {i}, starting with N0 = {N0}')  # debug
            S = (pc.R2S_lookup(R)  # S correction with nutation angle
                 - (Z - pc.R2Z_lookup(R)) * np.sin(np.radians(N0)))
            N1 = pc.S2N_lookup(S)  # update nutation with better S
            diff = N1 - N0
            # print(f'diff = {diff}')  # debug
            if np.all(np.abs(diff) < tol):  # deg
                break  # now we have S that satisfies the equality sufficiently
            else:
                N0 = N1
        # use the delta Z equation to calculate T since dZ is bigger than dR
        T = (Z - pc.S2Z_lookup(S)) / np.cos(np.radians(pc.S2N_lookup(S)))
        return np.vstack([Q, S, T])  # 3 x N array

    # %% QS transforms in 2D as projected from QST onto the focal surface

    @staticmethod
    def obsXY_to_QS(obsXY, cast=False):
        """
        On-shell condition is assumed, so T = 0
        input:  2 x N array
        output: 2 x N array, each column vector is QS
        """
        if cast:
            obsXY = typecast(obsXY)
        X, Y = obsXY[0, :], obsXY[1, :]
        Q = np.degrees(np.arctan2(Y, X))  # Y over X
        R = np.sqrt(np.square(X) + np.square(Y))
        S = pc.R2S_lookup(R)
        return np.vstack([Q, S])

    @staticmethod
    def obsXYZ_to_QS(obsXYZ, cast=False):
        '''obsXYZ -> QST -> QS, point can be off-shell and have nonzero T'''
        if cast:
            obsXYZ = typecast(obsXYZ)
        return PetalTransforms.obsXYZ_to_QST(obsXYZ)[:2, :]  # 2 x N array

    @staticmethod
    def QS_to_obsXYZ(QS, cast=False):
        """QS -> QST -> obsXYZ
        On-shell condition is assumed, so T = 0
        INPUT:  2 x N array, each column vector is QS
        OUTPUT: 3 x N array, each column vector is obsXYZ
        """
        if cast:
            QS = typecast(QS)
        QST = np.vstack([QS, np.zeros(QS.shape[1])])  # add zeros as T data
        return PetalTransforms.QST_to_obsXYZ(QST)

    @staticmethod
    def QS_to_obsXY(QS, cast=False):
        """QS -> QST -> obsXY
        On-shell condition is assumed, so T = 0
        INPUT:  2 x N array, each column vector is QS
        OUTPUT: 2 x N array, each column vector is obsXY
        """
        if cast:
            QS = typecast(QS)
        return PetalTransforms.QS_to_obsXYZ(QS)[:2, :]  # 2 x N array

    # %% flatXY transforms within petal local coordinates
    @staticmethod
    def ptlXYZ_to_flatXY(ptlXYZ, cast=False):
        """ptlXYZ -> QST -> QS -> flatXY
        INPUT:  3 x N or 2 x N array, each column vector is ptlXYZ or ptlXY
                in petal-local CS
        OUTPUT: 2 x N array, each column vector is flatXY, petal-local CS
        useful for seeding xy offsets on flattened focal surface
        """
        if cast:
            ptlXYZ = typecast(ptlXYZ)
        QS = PetalTransforms.obsXYZ_to_QS(ptlXYZ)  # assume location 3
        return PetalTransforms._QS_to_flatXY(QS)

    @staticmethod
    def flatXY_to_ptlXYZ(flatXY, cast=False):
        """flatXY -> QS -> QST -> ptlXYZ
        INPUT:  2 x N array, each column vector is flatXY, petal-local CS
        OUTPUT: 3 x N array, each column vector is ptlXYZ
                in petal-local CS
        """
        if cast:
            flatXY = typecast(flatXY)
        QS = PetalTransforms._flatXY_to_QS(flatXY)  # assume location 3
        return PetalTransforms.QS_to_obsXYZ(QS)

    @staticmethod
    def _QS_to_flatXY(QS, cast=False):
        """On-shell condition is assumed, global to global
        INPUT:  2 x N array, each column vector is QS
        OUTPUT: 2 x N array, each column vector is flatXY in CS5, reserved
        """
        if cast:
            QS = typecast(QS)
        Q_rad, S = np.radians(QS[0, :]), QS[1, :]
        X, Y = S * np.cos(Q_rad), S * np.sin(Q_rad)  # use S on shell as radius
        return np.vstack([X, Y])

    @staticmethod
    def _flatXY_to_QS(flatXY, cast=False):
        """On-shell condition is assumed, global to global
        INPUT:  2 x N array, each column vector is flatXY in CS5, reserved
        OUTPUT: 2 x N array, each column vector is QS
        """
        if cast:
            flatXY = typecast(flatXY)
        X, Y = flatXY[0, :], flatXY[1, :]
        Q = np.degrees(np.arctan2(Y, X))  # Y over X
        S = np.sqrt(np.square(X) + np.square(Y))  # R in flatXY is S on shell
        return np.vstack([Q, S])

    # %% petal transformations with 6 dof
    def ptlXYZ_to_obsXYZ(self, ptlXYZ, cast=False):
        """
        INPUT:  3 x N array, each column vector is ptlXYZ
        OUTPUT: 3 x N array, each column vector is obsXYZ
        """
        if cast:
            ptlXYZ = typecast(ptlXYZ)
        return self.R @ ptlXYZ + self.T  # forward transformation

    def obsXYZ_to_ptlXYZ(self, obsXYZ, cast=False):
        """
        INPUT:  3 x N array, each column vector is obsXYZ
        OUTPUT: 3 x N array, each column vector is ptlXYZ
        """
        if cast:
            obsXYZ = typecast(obsXYZ)
        return self.R.T @ (obsXYZ - self.T)  # backward transformation

    def QS_to_flatXY(self, QS, cast=False):
        """On-shell condition is assumed, global to local
        INPUT:  2 x N array, each column vector is QS, global
        OUTPUT: 2 x N array, each column vector is flatXY, petal-local
        """
        if cast:
            QS = typecast(QS)
        ptlXYZ = self.QS_to_ptlXYZ(QS)
        return self.ptlXYZ_to_flatXY(ptlXYZ)

    def flatXY_to_QS(self, flatXY, cast=False):
        """On-shell condition is assumed, local to global
        INPUT:  2 x N array, each column vector is flatXY, petal-local
        OUTPUT: 2 x N array, each column vector is QS, global
        """
        if cast:
            flatXY = typecast(flatXY)
        ptlXYZ = self.flatXY_to_ptlXYZ(flatXY)
        return self.ptlXYZ_to_QS(ptlXYZ)

    def flatXY_to_obsXY(self, flatXY, cast=False):
        """On-shell condition is assumed, local to global
        INPUT:  2 x N array, each column vector is flatXY, petal-local
        OUTPUT: 2 x N array, each column vector is obsXY, global
        """
        if cast:
            flatXY = typecast(flatXY)
        ptlXYZ = self.flatXY_to_ptlXYZ(flatXY)
        return self.ptlXYZ_to_obsXY(ptlXYZ)

    def obsXY_to_flatXY(self, obsXY, cast=False):
        """On-shell condition is assumed, global to local
        INPUT:  2 x N array, each column vector is obsXY, global
        OUTPUT: 2 x N array, each column vector is flatXY, petal-local
        """
        if cast:
            obsXY = typecast(obsXY)
        ptlXYZ = self.obsXY_to_ptlXYZ(obsXY)
        return self.ptlXYZ_to_flatXY(ptlXYZ)

    # %% composite transformations for convenience
    def ptlXYZ_to_QST(self, ptlXYZ, cast=False):
        """ptlXYZ -> obsXYZ -> QST -> QS
        INPUT:  3 x N array, each column vector is ptlXYZ, petal-local
        OUTPUT: 3 x N array, each column vector is QST
        """
        if cast:
            ptlXYZ = typecast(ptlXYZ)
        obsXYZ = self.ptlXYZ_to_obsXYZ(ptlXYZ)
        return self.obsXYZ_to_QST(obsXYZ)

    def QST_to_ptlXYZ(self, QST, cast=False):
        """QS -> QST -> obsXYZ -> ptlXYZ
        INPUT:  3 x N array, each column vector is QST
        OUTPUT: 3 x N array, each column vector is ptlXYZ, petal-local
        """
        if cast:
            QST = typecast(QST)
        obsXYZ = self.QST_to_obsXYZ(QST)
        return self.obsXYZ_to_ptlXYZ(obsXYZ)

    def ptlXYZ_to_QS(self, ptlXYZ, cast=False):
        """ptlXYZ -> obsXYZ -> QST -> QS
        INPUT:  3 x N array, each column vector is ptlXYZ, petal-local
        OUTPUT: 2 x N array, each column vector is QS
        """
        if cast:
            ptlXYZ = typecast(ptlXYZ)
        obsXYZ = self.ptlXYZ_to_obsXYZ(ptlXYZ)
        return self.obsXYZ_to_QS(obsXYZ)

    def QS_to_ptlXYZ(self, QS, cast=False):
        """QS -> QST -> obsXYZ -> ptlXYZ
        INPUT:  2 x N array, each column vector is QS
        OUTPUT: 3 x N array, each column vector is ptlXYZ, petal-local
        """
        if cast:
            QS = typecast(QS)
        obsXYZ = self.QS_to_obsXYZ(QS)
        return self.obsXYZ_to_ptlXYZ(obsXYZ)

    def ptlXYZ_to_obsXY(self, ptlXYZ, cast=False):
        """
        INPUT:  3 x N array, each column vector is ptlXYZ
        OUTPUT: 2 x N array, each column vector is obsXY
        """
        if cast:
            ptlXYZ = typecast(ptlXYZ)
        return self.ptlXYZ_to_obsXYZ(ptlXYZ)[:2, :]

    def obsXY_to_ptlXYZ(self, obsXY, cast=False):
        """
        INPUT:  2 x N array, each column vector is obsXY, global
        OUTPUT: 3 x N array, each column vector is ptlXYZ, local
        """
        if cast:
            obsXY = typecast(obsXY)
        QS = self.obsXY_to_QS(obsXY)
        return self.QS_to_ptlXYZ(QS)

# %% written but discouraged transformations

    # def QS_to_obsXY(self, QS, curved=None, cast=False):
    #     """
    #     input:  2 x N array, each column vector is QS
    #     output: 2 x N array, if curved=False, returns flatXY on focal surface
    #     """
    #     if cast:
    #         QS = typecast(QS)
    #     if curved is None:  # allow forcing flatXY despite curved instance
    #         curved = self.curved
    #     Q_rad, S = np.radians(QS[0, :]), QS[1, :]
    #     R = pc.S2R_lookup(S) if curved else S
    #     X = R * np.cos(Q_rad)
    #     Y = R * np.sin(Q_rad)
    #     return np.array([X, Y])

    # def ptlXY_to_obsXY(self, ptlXY, cast=False):
    #     """
    #     INPUT:  2 x N array, each column vector is ptlXY
    #     OUTPUT: 2 x N array, each column vector is obsXY
    #     On-shell condition is enforced because Z coordinates necessrily mix
    #     with XY when doing a (more accurate in principle) 3D transformation
    #     """
    #     if cast:
    #         ptlXY = typecast(ptlXY)
    #     R = np.sqrt(np.square(ptlXY[0, :]) + np.square(ptlXY[1, :]))
    #     ptlXYZ = np.vstack([ptlXY, pc.R2Z_lookup(R)])
    #     return self.ptlXYZ_to_obsXYZ(ptlXYZ)[:2, :]  # only take first 2 rows

    # def obsXY_to_ptlXY(self, obsXY, cast=False):
    #     """
    #     INPUT:  2 x N array, each column vector is obsXY
    #     OUTPUT: 2 x N array, each column vector is ptlXY
    #     On-shell condition is enforced
    #     """
    #     if cast:
    #         obsXY = typecast(obsXY)
    #     R = np.sqrt(np.square(obsXY[0, :]) + np.square(obsXY[1, :]))
    #     obsXYZ = np.vstack([obsXY, pc.R2Z_lookup(R)])
    #     return self.obsXYZ_to_ptlXYZ(obsXYZ)[:2, :]  # only take first 2 rows

    # def ptlXY_to_QS(self, ptlXY, cast=False):
    #     if cast:
    #         ptlXY = typecast(ptlXY)
    #     obsXY = self.ptlXY_to_obsXY(ptlXY)
    #     return self.obsXY_to_QS(obsXY)

    # def QS_to_ptlXY(self, QS, cast=False):
    #     if cast:
    #         QS = typecast(QS)
    #     obsXY = self.QS_to_obsXY(QS)
    #     return self.obsXY_to_ptlXY(obsXY)

    # def obsXY_to_flatXY(self, obsXY, cast=False):
    #     """Composite transformation, performs obsXY --> QS --> flatXY"""
    #     if cast:
    #         obsXY = typecast(obsXY)
    #     QS = self.obsXY_to_QS(obsXY)
    #     return self.QS_to_flatXY(QS)

    # def flatXY_to_obsXY(self, flatXY, cast=False):
    #     """Composite transformation, performs flatXY --> QS --> obsXY"""
    #     if cast:
    #         flatXY = typecast(flatXY)
    #     QS = self.flatXY_to_QS(flatXY)
    #     return self.QS_to_obsXY(QS)


# %% main
if __name__ == '__main__':
    """
    Unit test, choose location 4, which is 36 degrees of rotation ccw
    Taken from petal alignment data， PTL04
    T = [0.01261281, 0.068910657, 0.017850711]
    angles = [0.001382631,	-0.002945219,	35.9887675]
    dtb0_543 = [11.671,	23.542,	-81.893] in CS5
    [23.233778, 12.1450584, -81.9095268] in petal CS

    """
    trans = PetalTransforms(Tx=0.01261281, Ty=0.068910657, Tz=0.017850711,
                            alpha=0.001382631,
                            beta=-0.002945219,
                            gamma=0/180*np.pi)
    # obsXYZ_actual = np.array([11.671,	23.542,	-81.893]).reshape(3, 1)
    # ptlXYZ_actual = np.array([23.2337, 12.1450, -81.9095]).reshape(3, 1)
    # print(f'obsXYZ, actual: {obsXYZ_actual.T}\ntransformed: {obsXYZ.T}\n'
    #       f'ptlXYZ, actual: {ptlXYZ_actual.T}\ntransformed: {ptlXYZ.T}')
    trans = PetalTransforms(gamma=0/180*np.pi)
    # location 408
    ptlXYZ = np.array([346.797988, 194.710169, -17.737838]).reshape(3, 1)
    QST = np.array([29.312088, 398.257190, -0.2]).reshape(3, 1)  # 200 microns
    obsXYZ = trans.QST_to_obsXYZ(QST)
    print(f'obsXYZ = {obsXYZ.T}')
    QST = trans.obsXYZ_to_QST(obsXYZ)
    print(f'QST = {QST.T}')
    QS = trans.obsXYZ_to_QS(obsXYZ)
    print(f'QS = {QS.T}')
    obsXYZ = trans.QS_to_obsXYZ(QS)
    print(f'obsXYZ = {obsXYZ.T}')
    flatXY = trans.QS_to_flatXY(QS)
    print(f'flatXY = {flatXY.T}')
    QS = trans.flatXY_to_QS(flatXY)
    print(f'QS = {QS.T}')
    flatXY = trans.ptlXYZ_to_flatXY(ptlXYZ)
    print(f'flatXY = {flatXY.T}')
    ptlXYZ = trans.flatXY_to_ptlXYZ(flatXY)
    print(f'ptlXYZ = {ptlXYZ.T}')
    obsXYZ = trans.ptlXYZ_to_obsXYZ(ptlXYZ)
    print(f'obsXYZ = {obsXYZ.T}')
    ptlXYZ = trans.obsXYZ_to_ptlXYZ(obsXYZ)
    print(f'ptlXYZ = {ptlXYZ.T}')
    QS = trans.ptlXYZ_to_QS(ptlXYZ)
    print(f'QS = {QS.T}')
    ptlXYZ = trans.QS_to_ptlXYZ(QS)
    print(f'ptlXYZ = {ptlXYZ.T}')

    # new KEEPOUT_PTL definition
    ptlXY = np.array([[ 20.260, 410.189, 418.00,  406.3, 399.0, 384.0, 325.837, 20.260,  20.26, 420.0, 420.0, 20.26],
                      [  0.000,   0.000,  32.25 ,  89.0, 125.0, 167.5, 235.993, 13.978, 250.00 ,250.0 , -5.0, -5.00]])
    R = np.sqrt(np.square(ptlXY[0, :]) + np.square(ptlXY[1, :]))  # get Z
    Z = pc.R2Z_lookup(R)
    ptlXYZ = np.vstack([ptlXY, Z])
    print(f'flatXY KEEPOUT_PTL:\n{trans.ptlXYZ_to_flatXY(ptlXYZ)}')
    '''[[ 20.26005608 410.78494559 418.63957706 406.91006116 399.60662524
      384.58654953 326.28936763  20.26008312  20.27039953 361.45760424
      420.61560153  20.26005917]
     [  0.           0.          32.29934536  89.13363387 125.1900455
      167.75585168 236.32063497  13.97805735 250.12832589 215.15333585
       -5.00732859  -5.0000146 ]]'''

    # new KEEPOUT_GFA definition
    ptlXY = np.array([[295.569, 301.644, 303.538, 305.444, 307.547, 309.204, 320.770, 353.527],
                      [207.451, 201.634, 201.588, 201.296, 201.131, 199.968, 184.033, 207.831]])
    R = np.sqrt(np.square(ptlXY[0, :]) + np.square(ptlXY[1, :]))  # get Z
    Z = pc.R2Z_lookup(R)
    ptlXYZ = np.vstack([ptlXY, Z])
    print(f'flatXY KEEPOUT_GFA:\n{trans.ptlXYZ_to_flatXY(ptlXYZ)}')
    '''[[295.8911707  301.97618459 303.87536557 305.78637321 307.8951546
        309.55559758 321.1381771  354.04033488]
        [207.67712193 201.85604887 201.81205383 201.52163337 201.35868775
         200.19538472 184.24423152 208.13277865]]'''

    # additional QST tests
    obsXYZ = np.array([295.3328802, 207.2633873, -14.737207]).reshape(3, 1)
    print(f'Fitted obsXYZ = {obsXYZ.T}')
    QST = trans.obsXYZ_to_QST(obsXYZ)
    print(f'QST = {QST.T}')
    obsXYZ = trans.QST_to_obsXYZ(QST)
    print(f'obsXYZ = {obsXYZ.T}')

    # device_loc 541, fid
    obsXYZ = np.array([295.332581765, 207.26336747500002, -14.743447455]).reshape(3, 1)
    print(f'CMM obsXYZ = {obsXYZ.T}')
    QST = trans.obsXYZ_to_QST(obsXYZ)
    print(f'QST = {QST.T}')
    obsXYZ = trans.QST_to_obsXYZ(QST)
    print(f'obsXYZ = {obsXYZ.T}')
    
    # QS <-> obsXY transformation tests
    QS = np.array([65.55674, 411.2830]).reshape(2, 1)
    obsXY = trans.QS_to_obsXY(QS)
    print(f'obsXY = {obsXY.T}')
    flatXY = trans._QS_to_flatXY(QS)
    print(f'flatXY = {flatXY.T}')
    
    # QS obsXY consistency check
    ptlXYZ = np.array([346.797988, 194.710169, -17.737838]).reshape(3, 1)
    QS = trans.ptlXYZ_to_QS(ptlXYZ)
    print(f'QS = {QS.T}')
    obsXY1 = trans.ptlXYZ_to_obsXY(ptlXYZ)
    print(f'obsXY1 = {obsXY1.T}')
    obsXY2 = trans.QS_to_obsXY(QS)
    print(f'obsXY2 = {obsXY2.T}')