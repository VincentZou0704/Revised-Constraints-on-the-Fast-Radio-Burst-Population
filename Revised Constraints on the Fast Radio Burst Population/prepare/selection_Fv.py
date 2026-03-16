import FRBpopulation.prepare.CHIME_Lin as CL
from FRBpopulation.setup import *
import FRBpopulation.prepare.Fluence as fce
import time

def Ratio(lgfv, lgfv_max_th = lgfv_max_th0, n = n0):
    # grey zone difination
    eta_det = (lgfv - lgfv_min)/(lgfv_max_th - lgfv_min)
    eta = np.clip(eta_det, 0, 1)
    R = eta ** n
    return R


if __name__ == '__main__':

    start = time.time()
    plt.plot(lgfv_space, Ratio(lgfv_space))
    plt.show()
    end = time.time()
    print(end-start)




