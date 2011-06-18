import numpy
import scipy.special
from matplotlib import pyplot
from matplotlib import patches

def uniformOnSphere(r, n, m):
    if n > 1:
        z = numpy.random.randn(n, m)
        norms = (z**2).sum(axis=1)**0.5
        u = numpy.random.rand(n)**(1.0 / m)
        z /= norms[:, numpy.newaxis]
        z *= u[:, numpy.newaxis]
    else:
        z = numpy.random.randn(m)
        norms = (z**2).sum()**0.5
        u = numpy.random.rand()**(1.0 / m)
        z /= norms
        z *= u
    z *= r
    return z

def volume(r, m):
    return numpy.pi**(0.5*m) * r**m / scipy.special.gamma(0.5*m + 1)

def integrate1(mu, s, r, n):
    m = s.shape[0]
    z = numpy.random.randn(n, m)
    z /= s[numpy.newaxis,:]
    z += mu[numpy.newaxis,:]
    w = ((z**2).sum(axis=1) <= r**2).astype(float)
    g = w.sum() / n
    v = ((w - g)**2).sum() / (n * (n - 1))
    return g, v, w, z

def integrate2(mu, s, r, n):
    z = uniformOnSphere(r, n, s.shape[0])
    x = z - mu[numpy.newaxis, :]
    x *= s[numpy.newaxis, :]
    w = numpy.exp(-0.5 * numpy.sum(x**2, axis=1))
    w *= volume(r, s.shape[0]) * numpy.multiply.reduce(s) / (2.0 * numpy.pi)**(0.5 * s.shape[0])
    g = w.sum() / n
    v = ((w - g)**2).sum() / (n * (n - 1))
    return g, v, w, z

def integrate3(mu, s, r, n):
    g1, v1, w1, z1 = integrate1(mu, s, r, n)
    mask = (w1 == 0.0)
    n2 = mask.sum()
    if n2 < 0.50 * n:
        nTotal = n
        n1 = n - n2
        while n1 < n:
            ga, va, wa, za = integrate1(mu, s, r, n2)
            nTotal += n2
            w1[mask] = wa
            z1[mask] = za
            mask = (w1 == 0.0)
            n2 = mask.sum()
            n1 = n - n2
        g1 = w1.sum() / nTotal
        v1 = ((w1 - g1)**2).sum() / (nTotal * (nTotal - 1))
        return g1, v1, w1, z1
    g2, v2, w2, z2 = integrate2(mu, s, r, n2)
    if n2 < n - 1:
        v = (1.0 / (1.0 / v1 + 1.0 / v2))
        f1 = (1.0 / v1) * v
        f2 = (1.0 / v2) * v
    else:
        v = v2
        f1 = 0.0
        f2 = 1.0
    if not numpy.isfinite(v):
        print "v not finite:"
        print n-n2, f1, v1, g1
        print n2, f2, v2, g2
    w = w1 * f1
    w[mask] = w2 * f2
    z = z1.copy()
    z[mask] = z2
    g = g1*f1 + g2*f2
    return g, v, w, z

def scatter(mu, s, r, n=200):
    pyplot.figure(figsize=(18, 6))
    g1, v1, w1, z1 = integrate1(mu, s, r, n)
    g2, v2, w2, z2 = integrate2(mu, s, r, n)
    g3, v3, w3, z3 = integrate3(mu, s, r, n)
    def doPlot(ax, w, z):
        pyplot.scatter(z[:,0], z[:,1], c=w, alpha=0.1, linewidth=0)
        pyplot.axhline(mu[1], color='c')
        pyplot.axvline(mu[0], color='c')
        ax.add_patch(patches.Circle((0.0, 0.0), r, fill=False, linewidth=1, edgecolor="k"))
        ax.add_patch(patches.Ellipse((mu[0], mu[1]), 1.0/s[0], 1.0/s[1], fill=False, 
                                     linewidth=1, edgecolor="c"))
        pyplot.xlim(-5, 5)
        pyplot.ylim(-5, 5)
    doPlot(pyplot.subplot(1,3,1), w1, z1)
    doPlot(pyplot.subplot(1,3,2), w2, z2)
    doPlot(pyplot.subplot(1,3,3), w3, z3)
    print "%f +/- %f" % (g1, v1**0.5)
    print "%f +/- %f" % (g2, v2**0.5)
    print "%f +/- %f" % (g3, v3**0.5)

def sums(mu, s, r, k=50):
    pyplot.figure(figsize=(8, 6))
    x = numpy.arange(50, 1000, 50)
    y1 = numpy.zeros(x.shape, dtype=float)
    y1e = numpy.zeros(x.shape, dtype=float)
    y2 = numpy.zeros(x.shape, dtype=float)
    y2e = numpy.zeros(x.shape, dtype=float)
    y3 = numpy.zeros(x.shape, dtype=float)
    y3e = numpy.zeros(x.shape, dtype=float)
    for i, n in enumerate(x):
        for j in range(k):
            g1, v1, w1, z1 = integrate1(mu, s, r, n)
            g2, v2, w2, z2 = integrate2(mu, s, r, n)
            g3, v3, w3, z3 = integrate3(mu, s, r, n)
            y1[i] += g1
            y1e[i] += v1**0.5
            y2[i] += g2
            y2e[i] += v2**0.5
            y3[i] += g3
            y3e[i] += v3**0.5
    y1 /= k
    y2 /= k
    y3 /= k
    y1e /= k
    y2e /= k
    y3e /= k
    pyplot.fill_between(x, y1-y1e, y1+y1e, alpha=0.2, color='b')
    pyplot.fill_between(x, y2-y2e, y2+y2e, alpha=0.2, color='g')
    pyplot.fill_between(x, y3-y3e, y3+y3e, alpha=0.2, color='r')
    pyplot.plot(x, y1, "b", label="integrate1")
    pyplot.plot(x, y2, "g", label="integrate2")
    pyplot.plot(x, y3, "r", label="integrate3")
    pyplot.legend()

def comparison(m, n, k1=1000, k2=10):
    pyplot.figure()
    g = numpy.zeros(k1, dtype=float)
    vt1 = numpy.zeros(k1, dtype=float)
    vt2 = numpy.zeros(k1, dtype=float)
    mu = numpy.random.randn(k1, m)*1.0
    s = numpy.random.rand(k1, m)*1.0 + 1.0
    r = numpy.random.rand(k1)*3.0 + 1.0
    for i1 in xrange(k1):
        for i2 in xrange(k2):
            g1, v1, w1, z1 = integrate1(mu[i1], s[i1], r[i1], n)
            g2, v2, w2, z2 = integrate2(mu[i1], s[i1], r[i1], n)
            g[i1] += (g1 / v1 + g2 / v2) / (1.0 / v1 + 1.0 / v2)
            vt1[i1] += v1
            vt2[i1] += v2
    g /= k2
    vt1 /= k2
    vt2 /= k2
    rs = r * numpy.multiply.reduce(s, axis=1)**(1.0 / m)
    #pyplot.loglog(g, vt1, 'k,', alpha=0.5)
    pyplot.scatter(g, numpy.log10(vt1/vt2), c=rs, marker='d', alpha=0.5, linewidth=0)
    pyplot.axhline(0.0, color='k')
    pyplot.colorbar()
    pyplot.title("d=%d" % m)
    pyplot.ylabel("log10(Var(i1)/Var(i2))")
    pyplot.xlabel("integral")

if __name__ == "__main__":
    m = 2
    mu = numpy.random.randn(m)*1.5
    s = (numpy.random.rand(m)*2.0 + 1.0)
    s.sort()
    r = 2.0

    scatter(mu, s, r, 5000)
    sums(mu, s, r)
    #comparison(6, 100)
    pyplot.show()
