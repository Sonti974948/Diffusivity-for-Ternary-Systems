
import numpy as np
import math
from scipy.optimize import fsolve


#Acetic Acid Data 

p=101325 #Pa
Tb=391.2 #K
R=8.3144621
Tc=593 #K
Pc=57.87 #bar
w=0.447 
a=0.45724*((R*Tc)**2)/(Pc*100000)
b=0.0778*R*Tc/(Pc*100000)
Tr=Tb/Tc
alpha=(1+(0.37464+1.54226*w-0.26992*w*w)*(1-(Tr**0.5)))**2


guess=R*Tb/p

def PengR(v):    
    eqn=p-(R*Tb/(v-b))+(a*alpha/(v*(v+b)+b*(v-b)))    
    return eqn

Va=fsolve(PengR,guess) #m3/mol



def get_Da_w(Texp,musv,alsv,Msv):
    als=1.62
    Da_w=(7.4e-8)*Texp*np.sqrt(alsv*Msv/1000)/(musv/1000)/np.power((Va*(1e6)),0.6) #cm2/s

    return Da_w


def get_Da_org(Texp,musv,alsv,Msv):
    als=1.62
    Da_org=(7.4e-8)*Texp*np.sqrt(alsv*Msv/1000)/(musv/1000)/np.power(als*(Va*(1e6)),0.6) #cm2/s

    return Da_org





