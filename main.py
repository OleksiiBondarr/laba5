import matplotlib.mlab as ml
import pylab
import ll
import mag

print '6 proizvodnaya ', mag.df1(), '\n'
print 'supremum ', ll.sup(), '\n'
ll.newtonup()
print '\n'
ll.newtondown()
print '\n'
ll.lagr()
print '\n'
ll.splain()
print '\n'
ll.magoranta()
print '\n'
ll.gr()
'''
xmin = -20.0
xmax = 20.0
dx = 0.01
xlist = ml.frange(xmin, xmax, dx)
ylist = [f(x) for x in xlist]
pylab.plot(xlist, ylist)

pylab.show()
'''