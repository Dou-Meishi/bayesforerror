# -*- coding: utf-8 -*-
# generate some representative examples

import numpy as np
import matplotlib.pyplot as plt

import prde, dmsPDF             # personal module

def var_k():
    fn_fmt = dmsPDF.PDF.dirpath + '/' + "k-{k}_Q-{Q}"
    
    Q = .511
    c0ck = [1., 0., .11, 1.44, .25, -.41]
    abcent_id = [0, 1]

    def _k_ax(ax, k):
        pdfA = dmsPDF.PrDelta('A', Q, c0ck[:k+1], abcent_id)
        pdfB = dmsPDF.PrDelta('B', Q, c0ck[:k+1], abcent_id)
        pdfC = dmsPDF.PrDelta('C', Q, c0ck[:k+1], abcent_id)
        
        deltas = np.linspace(0, 2*pdfA.Q**(pdfA.k+1), 100)
        ax.plot(deltas, [pdfA(delta) for delta in deltas],
                '-', label="Set {}".format(pdfA.SET))
        ax.plot(deltas, [pdfB(delta) for delta in deltas],
                '--', label="Set {}".format(pdfB.SET))
        ax.plot(deltas, [pdfC(delta) for delta in deltas],
                ':', label="Set {}".format(pdfC.SET))
        ax.legend()
        return

    plt.rcParams.update({"font.size": 20, "lines.linewidth": 2})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 7.2))
    fig.suptitle("Q: {}".format(Q))
    
    k = 2
    ax1.set_title("k: {}".format(k))
    ax1.set_xlabel("$\Delta$")
    ax1.set_ylabel("pr($\Delta$)")
    _k_ax(ax1, k)

    k = 3
    ax2.set_title("k: {}".format(k))
    ax2.set_xlabel("$\Delta$")
    ax2.set_ylabel("pr($\Delta$)")
    _k_ax(ax2, k)

    fig.savefig(fn_fmt.format(k=23, Q=Q)+'.png')

    ax1.clear()
    ax2.clear()
    
    k = 4
    ax1.set_title("k: {}".format(k))
    ax1.set_xlabel("$\Delta$")
    ax1.set_ylabel("pr($\Delta$)")
    _k_ax(ax1, k)

    k = 5
    ax2.set_title("k: {}".format(k))
    ax2.set_xlabel("$\Delta$")
    ax2.set_ylabel("pr($\Delta$)")
    _k_ax(ax2, k)

    fig.savefig(fn_fmt.format(k=45, Q=Q)+'.png')
    return


def var_p():
    dataIX = {"68%":
              {'A': [.15, .097, .046, .022],
               'B': [.14, .091, .044, .022],
               'C': [.14, .096, .041, .019]},
              "95%":
              {'A': [.41, .20, .099, .041],
               'B': [.29, .17, .081, .040],
               'C': [.47, .26, .100, .043]}
              }

    plt.rcParams.update({"font.size": 20})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 7.2))
    fig.suptitle("Q: {}".format(.511))

    ax1.set_xlabel('k')
    ax1.set_ylabel("$d_k^p$")
    ax2.set_xlabel('k')
    ax2.set_ylabel("$d_k^p$")

    ks = [2, 3, 4, 5]
    ax1.set_xticks(ks)
    ax2.set_xticks(ks)

    p = "68%"
    ax1.set_title("p: {}".format(p))

    ax1.plot(ks, dataIX[p]['A'], 'X', label="Set: A")
    ax1.plot(ks, dataIX[p]['B'], 'P', label="Set: B")
    ax1.plot(ks, dataIX[p]['C'], 'd', label="Set: C")
    ax1.legend()
    
    p = "95%"
    ax2.set_title("p: {}".format(p))

    ax2.plot(ks, dataIX[p]['A'], 'X', label="Set: A")
    ax2.plot(ks, dataIX[p]['B'], 'P', label="Set: B")
    ax2.plot(ks, dataIX[p]['C'], 'd', label="Set: C")
    ax2.legend()

    fig.savefig(dmsPDF.PDF.dirpath+"/var_p.png")
    return


def main():
    var_k()
    # var_p()

    return


if __name__ == '__main__':
    import cProfile
    cProfile.run("main()", filename="first-report.profile")
