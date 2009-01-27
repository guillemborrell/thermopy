# -*- coding: utf-8 -*-
from __future__ import division
from units import Pressure,Temperature,Enthalpy
from numpy import array,sum,sqrt


class Water(object):
    """Taken from

    The International Association for the Properties of Water and
    Steam. Lucerne, Switzerland. August 2007. Revised Release on the
    IAPWS Industrial Formulation 1997 for the Thermodynamic Properties
    of Water and Steam.
    
    Functions implemented:
    
    Saturation line
    h(p,T) region1, region2 and no warnings yet
    """

    R=0.461526 #kJ/(kg K)
    Tc=647.096 #Triple point temperature (K)
    pc=22.064  #Triple point pressure (MPa)
    rc=322     #Triple point density kg/m3

    data={
        'n':array([0.11670521452767E4,   
                   - 0.72421316703206E6,
                   - 0.17073846940092E2,
                   0.12020824702470E5  ,
                   - 0.32325550322333E7,
                   0.14915108613530E2  ,
                   - 0.48232657361591E4,
                   0.40511340542057E6  ,
                   - 0.23855557567849  ,
                   0.65017534844798E3],'d'),
        'pT_datar1':array([
                [0 ,-2 , 0.14632971213167    ],
                [0 ,-1 ,-0.84548187169114    ],
                [0 , 0 ,-0.37563603672040E1  ],
                [0 , 1 , 0.33855169168385E1  ],
                [0 , 2 ,-0.95791963387872    ],
                [0 , 3 , 0.15772038513228    ],
                [0 , 4 ,-0.16616417199501E-1 ],
                [0 , 5 , 0.81214629983568E-3 ],
                [1 ,-9 , 0.28319080123804E-3 ],
                [1 ,-7 ,-0.60706301565874E-3 ],
                [1 ,-1 ,-0.18990068218419E-1 ],
                [1 , 0 ,-0.32529748770505E-1 ],
                [1 , 1 ,-0.21841717175414E-1 ],
                [1 , 3 ,-0.52838357969930E-4 ],
                [2 ,-3 ,-0.47184321073267E-3 ],
                [2 , 0 ,-0.30001780793026E-3 ],
                [2 , 1 , 0.47661393906987E-4 ],
                [2 , 3 ,-0.44141845330846E-5 ],
                [2 , 17,-0.72694996297594E-15],
                [3 ,-4 ,-0.31679644845054E-4 ],
                [3 , 0 ,-0.28270797985312E-5 ],
                [3 , 6 ,-0.85205128120103E-9 ],
                [4 ,-5 ,-0.22425281908000E-5 ],
                [4 ,-2 ,-0.65171222895601E-6 ],
                [4 , 10,-0.14341729937924E-12],
                [5 ,-8 ,-0.40516996860117E-6 ],
                [8 ,-11,-0.12734301741641E-8 ],
                [8 ,-6 ,-0.17424871230634E-9 ],
                [21,-29,-0.68762131295531E-18],
                [23,-31, 0.14478307828521E-19],
                [29,-38, 0.26335781662795E-22],
                [30,-39,-0.11947622640071E-22],
                [31,-40, 0.18228094581404E-23],
                [32,-41,-0.93537087292458E-25]],'d'),
        'pT_datar20':array([
                [0 , -0.96927686500217E1 ], 
                [1 , 0.10086655968018E2  ],
                [-5, -0.56087911283020E-2],
                [-4, 0.71452738081455E-1 ],
                [-3, -0.40710498223928   ],
                [-2, 0.14240819171444E1  ],
                [-1, -0.43839511319450E1 ],
                [2 , -0.28408632460772   ],
                [3 , 0.21268463753307E-1 ]],'d'),
        'pT_datar2r':array([
                [1 ,0 ,-0.17731742473213E-2],
                [1 ,1 ,-0.17834862292358E-1],
                [1 ,2 ,-0.45996013696365E-1],
                [1 ,3 ,-0.57581259083432E-1],
                [1 ,6 ,-0.50325278727930E-1],
                [2 ,1 ,-0.33032641670203E-4],
                [2 ,2 ,-0.18948987516315E-3],
                [2 ,4 ,-0.39392777243355E-2],
                [2 ,7 ,-0.43797295650573E-1],
                [2 ,36,-0.26674547914087E-4],
                [3 ,0 , 0.20481737692309E-7],
                [3 ,1 , 0.43870667284435E-6],
                [3 ,3 ,-0.32277677238570E-4],
                [3 ,6 ,-0.15033924542148E-2],
                [3 ,35,-0.40668253562649E-1],
                [4 ,1 ,-0.78847309559367E-9],
                [4 ,2 , 0.12790717852285E-7],
                [4 ,3 , 0.48225372718507E-6],
                [5 ,7 , 0.22922076337661E-5],
                [6 ,3 ,-0.16714766451061E-1],
                [6 ,16,-0.21171472321355E-2],
                [6 ,35,-0.23895741934104E2 ],
                [7 ,0 ,-0.59059564324270E-1],
                [7 ,11,-0.12621808899101E-5],
                [7 ,25,-0.38946842435739E-1],
                [8 ,8 , 0.11256211360459E-1],
                [8 ,36,-0.82311340897998E1 ],
                [9 ,13, 0.19809712802088E-7],
                [10,4 , 0.10406965210174E-1],
                [10,10,-0.10234747095929E-1],
                [10,14,-0.10018179379511E-8],
                [16,29,-0.80882908646985E-1],
                [16,50, 0.10693031879409   ],
                [18,57,-0.33662250574171   ],
                [20,20, 0.89185845355421E-2],
                [20,35, 0.30629316876232E-1],
                [20,48,-0.42002467698208E-5],
                [21,21,-0.59056029685639E-2],
                [22,53, 0.37826947613457E-5],
                [23,39,-0.12768608934681E-1],
                [24,26, 0.73087610595061E-2],
                [24,40, 0.55414715350778E-1],
                [24,58,-0.94369707241210E-6]],'d'),
        }
    def psat(self,T):
        """
        Returns the saturation pressure of water at a given temperature.
        
        Remember that temperature must be between 273.15K (triple point)
        and 647.096K (critical point)
        
        Temperatures in K, Pressures in Pa.
        
        >>> w = Water()
        >>> w.psat(300)
        3536.5894130130105
        >>> w.psat(130)
        Traceback (most recent call last):
          File "/usr/lib/python2.5/doctest.py", line 1228, in __run
            compileflags, 1) in test.globs
          File "<doctest __main__.Water.psat[2]>", line 1, in <module>
            w.psat(130)
          File "iapws.py", line 153, in psat
            'No saturation pressure for this temperature')
        ValueError: No saturation pressure for this temperature
        >>> w.psat(700)
        Traceback (most recent call last):
          File "/usr/lib/python2.5/doctest.py", line 1228, in __run
            compileflags, 1) in test.globs
          File "<doctest __main__.Water.psat[3]>", line 1, in <module>
             w.psat(700)
          File "iapws.py", line 146, in psat
            'No saturation pressure for this temperature')
        ValueError: No saturation pressure for this temperature
        """
        if T < 273.15 or T > 647.096:
            raise ValueError(
                'No saturation pressure for this temperature')
        
        n=self.data['n']
        v=T+n[8]/(T-n[9])
        A=v**2+n[0]*v+n[1]
        B=n[2]*v**2+n[3]*v+n[4]
        C=n[5]*v**2+n[6]*v+n[7]
        return Pressure(((2*C)/(-B+sqrt(B**2-4*A*C)))**4).unit('MPa')


    def Tsat(self,p):
        """
        Returns the saturation temperature of water at a given pressure.
        
        Remember that pressure must be between 0.000611213MPa (triple
        point) and 22.064MPa (critical point)
        
        Temperatures in K, pressures in MPa.

        >>> w = Water()
        >>> w.Tsat(100000)
        372.75591861133762
        >>> w.Tsat(1200000)
        461.11464162130039
        >>> w.Tsat(100)
        Traceback (most recent call last):                                 
          File "/usr/lib/python2.5/doctest.py", line 1228, in __run        
            compileflags, 1) in test.globs                                 
          File "<doctest __main__.Water.Tsat[3]>", line 1, in <module>     
            w.Tsat(100)                                                    
          File "iapws.py", line 193, in Tsat                               
            raise ValueError('No saturation temperature for this pressure')
        ValueError: No saturation temperature for this pressure
        >>> w.Tsat(101325)
        373.12430000048056
        """
        
        p = Pressure(p).MPa

        if p < 0.000611213 or p > 22.064:
            raise ValueError('No saturation temperature for this pressure')

        n=self.data['n']
        beta=p**0.25
        E=beta**2+n[2]*beta+n[5]
        F=n[0]*beta**2+n[3]*beta+n[6]
        G=n[1]*beta**2+n[4]*beta+n[7]
        D=(2*G)/(-F-sqrt(F**2-4*E*G))
        
        return Temperature(0.5*(n[9]+D-sqrt((n[9]+D)**2-4*(n[8]+n[9]*D))))


    def h(self,p,T):
        """
        Returns specific enthalpy (J/kg) for a given pressure (Pa) and
        Temperature (K).
        
        >>> w = Water()
        >>> w.h(3000000,300)
        115331.27302143839
        >>> w.h(3500,300)
        2549911.4508400201

        There are also error codes        

        Results checked against the reference.
        """
        p = Pressure(p).MPa
        
        # region 1 implementation
        if p >= self.psat(T).MPa:
        # Liquid water (pressure over saturation pressure)
            pi=p/16.53
            tau=1386/T
    
            raw_data=self.data['pT_datar1']
            I=raw_data[:,0]
            J=raw_data[:,1]
            n=raw_data[:,2]

            return Enthalpy(
                self.R*T*tau*sum((n*(7.1-pi)**I)*J*((tau-1.222)**(J-1)))).unit('kJkg')

        if p < self.psat(T).MPa:
            # steam, pressure under saturation pressure
            pi=p
            tau=540/T 
            
            raw_data0=self.data['pT_datar20']
            J0=raw_data0[:,0]
            n0=raw_data0[:,1]
            
            raw_datar=self.data['pT_datar2r']
            I=raw_datar[:,0]
            J=raw_datar[:,1]
            n=raw_datar[:,2]
            
            g0_tau=sum(n0*J0*tau**(J0-1))
            gr_tau=sum(n*pi**I*J*(tau-0.5)**(J-1))
            
            return Enthalpy(self.R*T*tau*(g0_tau+gr_tau)).unit(
                'kJkg')


    def T_ph(self,p,h):
        """
        Returns the temperature (K) given the pressure (MPa) and
        specific enthalpy (kJ/kg).  Only region 2a implemented
        (p<4MPa) (Reimplement).

        >>> w = Water()
        >>> w.T_ph(3,500)
        (391.79850876242563, 4.1313215739117547e+21)
        >>> w.T_ph(3,4000)
        (-14923984.403552555, 1010.7757656610391)
        >>> w.T_ph(0.001,3000)
        (-103213.84623418792, 534.43324138173193)
        >>> w.T_ph(3,500)
        (391.79850876242563, 4.1313215739117547e+21)
        >>> w.T_ph(80,500)
        (378.10862587940176, -6.0291236598288914e+28)
        >>> w.T_ph(80,1500)
        (611.04122940266461, -5.5726211553407323e+22)
        """
        eta = h/2500.
        
        raw_data=array([
                [1 ,0,0 ,-0.23872489924521E3 ],
                [2 ,0,1 , 0.40421188637945E3 ],  
                [3 ,0,2 , 0.11349746881718E3 ],  
                [4 ,0,6 ,-0.58457616048039E1 ],
                [5 ,0,22,-0.15285482413140E-3],
                [6 ,0,32,-0.10866707695377E-5],
                [7 ,1,0 ,-0.13391744872602E2 ],
                [8 ,1,1 , 0.43211039183559E2 ],  
                [9 ,1,2 ,-0.54010067170506E2 ],
                [10,1,3 , 0.30535892203916E2 ],  
                [11,1,4 ,-0.65964749423638E1  ],
                [12,1,10, 0.93965400878363E-2 ],  
                [13,1,32, 0.11573647505340E-6 ],  
                [14,2,10,-0.25858641282073E-4 ],
                [15,2,32,-0.40644363084799E-8 ],
                [16,3,10, 0.66456186191635E-7 ],  
                [17,3,32, 0.80670734103027E-10],  
                [18,4,32,-0.93477771213947E-12],
                [19,5,32, 0.58265442020601E-14],  
                [20,6,32,-0.15020185953503E-16]],'d')
        
        i=raw_data[:,0]
        I=raw_data[:,1]
        J=raw_data[:,2]
        n=raw_data[:,3]
        
        raw_data2=array([
                [1 ,0,0 , 0.10898952318288E4 ],  
                [2 ,0,1 , 0.84951654495535E3 ],  
                [3 ,0,2 ,-0.10781748091826E3 ],
                [4 ,0,3 , 0.33153654801263E2 ],  
                [5 ,0,7 ,-0.74232016790248E1 ],
                [6 ,0,20, 0.11765048724356E2 ],  
                [7 ,1,0 , 0.18445749355790E1 ],  
                [8 ,1,1 ,-0.41792700549624E1 ],
                [9 ,1,2 , 0.62478196935812E1 ],  
                [10,1,3 ,-0.17344563108114E2 ],
                [11,1,7 ,-0.20058176862096E3 ],
                [12,1,9 , 0.27196065473796E3 ],  
                [13,1,11,-0.45511318285818E3 ],
                [14,1,18, 0.30919688604755E4 ],  
                [15,1,44, 0.25226640357872E6 ],  
                [16,2,0 ,-0.61707422868339E-2],
                [17,2,2 ,-0.31078046629583   ],       
                [18,2,7 , 0.11670873077107E2 ],  
                [19,2,36, 0.12812798404046E9 ],  
                [20,2,38,-0.98554909623276E9 ],
                [21,2,40, 0.28224546973002E10],  
                [22,2,42,-0.35948971410703E10],
                [23,2,44, 0.17227349913197E10],  
                [24,3,24,-0.13551334240775E5 ],
                [25,3,44, 0.12848734664650E8 ],  
                [26,4,12, 0.13865724283226E1 ],  
                [27,4,32, 0.23598832556514E6 ],  
                [28,4,44,-0.13105236545054E8 ],
                [29,5,32, 0.73999835474766E4 ],  
                [30,5,36,-0.55196697030060E6 ],
                [31,5,42, 0.37154085996233E7 ],  
                [32,6,34, 0.19127729239660E5 ],  
                [33,6,44,-0.41535164835634E6 ],
                [34,7,28,-0.62459855192507E2 ]],'d')
        
        eta2 = h/2000.
        i2=raw_data2[:,0]
        I2=raw_data2[:,1]
        J2=raw_data2[:,2]
        n2=raw_data2[:,3]
        
        return (sum(n*p**I*(eta+1)**J),sum(n2*p**I2*(eta2-2.1)**J2))
    

def test():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    test()

