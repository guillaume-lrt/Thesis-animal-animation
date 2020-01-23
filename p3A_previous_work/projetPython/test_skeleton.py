from skimage.morphology import skeletonize
from skimage.morphology import local_maxima
from matplotlib import colors
import matplotlib.pyplot as plt
from skimage.util import invert
from scipy import misc
import cv2
import numpy as np

def rgb_to_gray(img):
        R = np.array(img[:, :, 0])
        G = np.array(img[:, :, 1])
        B = np.array(img[:, :, 2])

        R = (R *.299)
        G = (G *.587)
        B = (B *.114)

        Avg = (R+G+B)

        grayImage = Avg
        print(grayImage[100:105,115:117])
        return grayImage

def plot(imag):
    # Invert the horse image
    # image = plt.imread("cheval.png", 'L')/255
    image = rgb_to_gray(1-plt.imread(imag))
    print(image.shape)
    print(image[100:105,115:117])
    print(type(image))

    # perform skeletonization
    # skeleton = skeletonize(image)
    skeleton = local_maxima(image)

    # display results
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4),
                             sharex=True, sharey=True,
                             subplot_kw={'adjustable': 'box'})

    ax = axes.ravel()

    ax[0].imshow(image, cmap=plt.cm.gray)
    ax[0].axis('off')
    ax[0].set_title('original', fontsize=20)

    ax[1].imshow(skeleton, cmap=plt.cm.gray)
    ax[1].axis('off')
    ax[1].set_title('skeleton', fontsize=20)

    fig.tight_layout()
    plt.show()

# plot("elephant.jpg")
# plot("elephant.png")
plot("elephant3.jpg")
plot("0_rgba.png")
plot("1_rgba.png")
plot("5_rgba.png")
plot("6_rgba.png")
