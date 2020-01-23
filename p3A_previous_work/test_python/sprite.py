#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 17:50:38 2018

@author: coraliedavid
"""
import pygame
from constantes import *
import cv2

class Sprite:
    """ Displays 3D objects on a Pygame screen """

    def __init__(self, width, height, v_n, v_t, names = []):
        self.width = width
        self.height = height
        self.images = []
        for name in names:
            #im = pygame.image.load(name).convert_alpha()
            #im.set_color_key(BLACK)
            self.images.append(im)
        
        self.vect_n = v_n
        self.vect_t = v_t
        self.pos = 0
        
    def display_image():
        M = np.zeros()
        M[0] = v_t
        M[1] = np.cross(v_n, v_t)
        dst = cv2.warpAffine(images[self.pos], M,(cols,rows))
        plt.subplot(121),plt.imshow(img),plt.title('Input')
        plt.subplot(122),plt.imshow(dst),plt.title('Output')
        plt.show()