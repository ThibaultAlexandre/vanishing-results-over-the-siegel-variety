"""``p``-small weights for ``Sp_{2g}``, twisted by ``-w_0``.

A dominant weight ``chi`` is ``p``-small when ``|<chi + rho, alpha>| <= p`` for
every positive root ``alpha`` of ``C_g``.
"""

from sage.all import Integer, QQ


def psmall(data, kmax, directory="save", write=True):
    """Return, and optionally save, the ``p``-small weights up to ``kmax``."""
    rho = sum(data.phi_G) / 2
    rho_vec = tuple(QQ(rho[j]) for j in range(data.g))
    phi_vec = [tuple(Integer(a[j]) for j in range(data.g)) for a in data.phi_G]

    res = []
    for char in data.buildWeights(0, kmax):
        ok = True
        for alpha in phi_vec:
            prod = sum((char[j] + rho_vec[j]) * alpha[j] for j in range(data.g))
            if abs(prod) > data.p:
                ok = False
                break
        if ok:
            res.append(data.changeConvention(char))

    if write:
        path = "%s/g%dp%d_psmall.txt" % (directory, data.g, data.p)
        with open(path, "w") as f:
            for char in res:
                f.write("".join("%s " % x for x in char) + "\n")
        print("p-small weights saved in " + path)
    return res
