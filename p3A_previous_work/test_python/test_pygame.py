import sys, pygame
from pygame.locals import *

pygame.init()

size = width, height = 800, 700
speed = [2, 2]
black = 0, 0, 0
white = 255, 255, 255

screen = pygame.display.set_mode(size)

ball = pygame.image.load("ball.jpg").convert_alpha()
    
ballrect = ball.get_rect()
perso_x, perso_y = 0, 0

continuer = 1

while continuer:
    for event in pygame.event.get():
        if event.type == QUIT: 
            
            continuer = 0
        if event.type == MOUSEBUTTONDOWN: 
            print(event.button)
            if event.button == 1:	#Si clic gauche
                perso_x = event.pos[0]
                perso_y = event.pos[1]

    screen.fill(white)
    screen.blit(ball, (perso_x, perso_y))
    pygame.display.flip()



sys.exit()