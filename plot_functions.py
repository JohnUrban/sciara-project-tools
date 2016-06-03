from matplotlib import pyplot as plt

def scatterplot(X, Y, vmin, vmax, size=False, names=False, saveas=False, do_nothing=False, color=None, alpha=0.5, grayscale=False, colorbar=False, cmap='seismic'):
    if not size:
        size = 10
    if color is None:
        color = 'b'
        colorbar=False ##ENSURE
        vmin=1
        vmax=1
    plt.scatter(x=X, y=Y, s=size, c=color, vmin=vmin, vmax=vmax, alpha=alpha)#, cmap=cmap) 
##    plt.scatter(x=X, y=Y, s=size, c=color, cmap=cmap) 
    if grayscale:
        plt.gray()
    if names:
        for i in range(len(NR)):
            plt.annotate(names[i], (X[i], Y[i]), size='xx-small')
    if colorbar:
        plt.colorbar()
    if saveas is not None and (saveas.endswith(".pdf") or saveas.endswith(".jpg")):
            plt.savefig(saveas)
    elif do_nothing:
        return
    else:
        plt.show()
