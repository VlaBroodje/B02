import math
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict


Cycle = 0

class MassEstimate:
    Maxmass = 9593.7
    RangeMaxmass = 5116.72
    Emptymass = 5852.16
    RangeEmptymass = 7192.62

while Cycle<100:

    deltaT = 0
    pISA = 101325
    rhoISA = 1.225
    tempISA = 288.15+deltaT

    Clmax = 1.6
    Cllanding = 2.4
    ClTakeoff = 2
    rhoLanding = 1.1646

    def rho(P, T):
        rho = P/(T*287)
        return rho

    class tlar:
        LengthTakeoff = 1250
        minclimbrate = 0.5
        MfracClimb = 0.98
        Mach = 0.68
        Hcruise = 12496.8
        MfracCcruise = 0.95
        LengthLanding = 620
        MfracLanding = 0.91
        MfracTakeoff = 0.995
        VAproach = 45.94
        Soundspeed = 295
        ClimbGradient119 = 3.2
        ClimbGradient121a = 0
        ClimbGradient121b = 2.4
        ClimbGradient121c = 1.2
        ClimbGradient121d = 2.1

    class DragPolar:
        AR = 9
        frictionCoefficient = 0.0035
        CD0 = 0.0175
        parasiteDrag = 0.0075
        spanEfficiency = 0.97
        Oswald = 1/(math.pi*AR*parasiteDrag+(1/spanEfficiency))
        Bypass = 5
        jeteff = ((tlar.Mach*tlar.Soundspeed)/((Bypass**-0.19)*22))*1/44

    
    class Wing:
        a = 2

    pCruise = 17873.87
    tempCruise = 216.65+deltaT
    rhoCruise = pCruise/(287*tempCruise)
    vCruise = tlar.Mach*tlar.Soundspeed
    pTakeoff = 95459.79
    rhoTakeoff = 1.116
    ClstallTakeoff = ((1/1.13)**2)*ClTakeoff






    class Point:
        def __init__(self, x, y):
            self.x = x
            self.y = y
        

    ###0/10/30    cruise/takeoff/landing
    def config(deflection, gear):
        CD0 = DragPolar.CD0
        Oswald = DragPolar.Oswald
        if gear == 1:
            CD0 += 0.02
        CD0 += 0.0013*deflection
        Oswald +=0.0046*deflection
        return CD0, Oswald

    def ClmaxGradient(Cd0, Oswald, AR):
        return math.sqrt(Cd0*math.pi*AR*Oswald)

    def lapse(p,t, M):
        tbreak = 1.06
        pT = p*((1+0.2*(M**2))**3.5)
        deltaT = pT/pISA
        tempT = t*(1+0.2*(M**2))
        thetatT = tempT/tempISA
        if(thetatT<tbreak):
            lapse = deltaT*(1-(0.43+0.014*DragPolar.Bypass)*math.sqrt(M))
        else:
            lapse = deltaT*(1-(0.43+0.014*DragPolar.Oswald)*math.sqrt(M))-(3*(thetatT-tbreak))/(1.5+M)
        return lapse

    def minSpeed(Mfrac, Vapproach, Cl):
        values = []
        i = 0
        while(i<35):
            Loading = (1/Mfrac)*(1.225/2)*((Vapproach/1.23)**2)*Cl
            values.append(Loading)
            i+=1
        return values

    def LandingLength(Mfrac, fieldlength, rho, Cl):
        values = []
        i = 0
        while(i<35):
            Loading = (1/Mfrac)*(fieldlength/0.45)*(rho*Cl/2)
            values.append(Loading)
            i+=1
        return values

    def CruiseSpeed (Mfrac, tlapse, Cd0, rhoCruise, VCruise, Aspect, Oswald):
        values = []
        i = 0
        L = 1200
        while(L<3500):
            Loading = (Mfrac/tlapse)*((Cd0*0.5*rhoCruise*(VCruise**2)) / (Mfrac*L) + (Mfrac*L) / (math.pi * Aspect * Oswald * 0.5 * rhoCruise * (VCruise ** 2)))
            L+=100
            #print(L)
            values.append(Point(L, Loading))
            i+=1
        return values

    def ClimbRate (Mfrac, climbRate, rhoClimb, Cd0, Aspect, Oswald):
        values = []
        i = 0
        L = 100
        Cl = math.sqrt(Cd0*math.pi*Aspect*Oswald)
        while(i<35):
            VClimb = math.sqrt(2 * L * Mfrac /(rhoClimb * Cl))
            Mach = VClimb / math.sqrt(1.4 * 287 * tempCruise)
            tlapse = lapse(pCruise, tempCruise, Mach)
            Loading = (Mfrac/tlapse)*((math.sqrt(((climbRate**2)*rhoClimb)/(2*Mfrac*L)*math.sqrt(Cd0*math.pi*Aspect*Oswald)))+2*math.sqrt(Cd0/(math.pi*Aspect*Oswald)))
            values.append(Point(L, Loading))
            L+=100
            i+=1
        return values

    def ClimbGradient (ClimbGradient, deflection, Aspect, gear):
        Cd0 = config(deflection, gear)[0]
        Oswald = config(deflection, gear)[1]
        rho = 1.225
        Mfrac = 1
        values = []
        i = 0
        L = 100
        Cl = math.sqrt(config(deflection,gear)[0] * math.pi * Aspect * config(30,1)[1])
        #print("PEDERAST",Cl)
        while (i < 35):
            VClimb = math.sqrt(2 * L * Mfrac / (rho * Cl))
            Mach = VClimb / math.sqrt(1.4 * 287 * tempISA)
            tlapse = lapse(pISA, tempISA, Mach)
            #print(tlapse)
            Loading = (Mfrac/tlapse)*((ClimbGradient/100)+2*math.sqrt(Cd0/(math.pi*Aspect*Oswald)))
            values.append(Point(L, Loading))
            L += 100
            i += 1
        return values

    def Takeoff (Cl, Oswald, Aspect, rho):
        values = []
        Mfrac = 1
        L = 100
        i = 0
        while (i<35):
            VClimb = math.sqrt(2 * L * Mfrac / (rho * Cl))
            Mach = VClimb / math.sqrt(1.4 * 287 * tempISA)
            tlapse = lapse(pISA, tempISA, Mach)
            Loading = (Mfrac/tlapse)*(1.15*math.sqrt(2*L/(1250*0.85*9.80665*math.pi*Aspect*Oswald))+(88)/1250)
            values.append(Point(L, Loading))
            i+=1
            L+=100
        return values

    Loading = []
    L = 100
    i = 0
    while(i<35):
        Loading.append(L)
        L+=100
        i+=1

    bruh = []
    b = 0
    i = 0
    while(i<35):
        bruh.append(b)
        b+=0.02
        i+=1

    def PointToNum(point):
        nums = []
        i = 0
        while(i<len(point)):
            nums.append(point[i].y)
            i+=1
        return nums


    def PointToNumx(point):
        nums = []
        i = 0
        while(i<len(point)):
            nums.append(point[i].x)
            i+=1
        return nums

    def searchXreturnY(list, x):
        found = False
        y = 0
        i = 0
        while (i<len(list) and found == False):
            if(list[i].x>x):
                x1 = list[i-1].x
                x2 = list[i].x
                found = True
                slope = (list[i].y-list[i-1].y)/(x2-x1)
                y = list[i-1].y+slope*(x-x1)
            i+=1

        return y



    cruiseSpeed = CruiseSpeed(tlar.MfracCcruise, lapse(pCruise, tempCruise, tlar.Mach), DragPolar.CD0, rhoCruise, vCruise, DragPolar.AR, DragPolar.Oswald)
    climbRate = ClimbRate(tlar.MfracClimb, tlar.minclimbrate, rhoCruise, DragPolar.CD0, DragPolar.AR, DragPolar.Oswald)
    Cs119 = ClimbGradient(tlar.ClimbGradient119, 30, DragPolar.AR, 1)
    Cs121a = ClimbGradient(tlar.ClimbGradient121a, 10, DragPolar.AR, 1)
    Cs121b = ClimbGradient(tlar.ClimbGradient121b, 10, DragPolar.AR, 0)
    Cs121c = ClimbGradient(tlar.ClimbGradient121c, 0, DragPolar.AR, 0)
    Cs121d = ClimbGradient(tlar.ClimbGradient121d, 30, DragPolar.AR, 0)
    TakeoffLength = Takeoff(ClstallTakeoff, DragPolar.Oswald, DragPolar.AR, rhoTakeoff)

    Points = []
    Points.append(cruiseSpeed)
    Points.append(climbRate)
    Points.append(Cs119)
    Points.append(Cs121a)
    Points.append(Cs121b)
    Points.append(Cs121c)
    Points.append(Cs121d)
    Points.append(TakeoffLength)



    MinSpeedy = minSpeed(tlar.MfracTakeoff, tlar.VAproach, Cllanding)
    landingLengthy = LandingLength(tlar.MfracLanding, tlar.LengthLanding, rhoLanding, Cllanding)
    cruiseSpeedy = PointToNum(cruiseSpeed)
    climbRatey = PointToNum(climbRate)
    Cs119y = PointToNum(Cs119)
    Cs121ay = PointToNum(Cs121a)
    Cs121by = PointToNum(Cs121b)
    Cs121cy = PointToNum(Cs121c)
    Cs121dy = PointToNum(Cs121d)
    TakeoffLengthy = PointToNum(TakeoffLength)


    cruiseSpeedx = PointToNumx(cruiseSpeed)
    climbRatex = PointToNumx(climbRate)
    Cs119x = PointToNumx(Cs119)
    Cs121ax = PointToNumx(Cs121a)
    Cs121bx = PointToNumx(Cs121b)
    Cs121cx = PointToNumx(Cs121c)
    Cs121dx = PointToNumx(Cs121d)
    TakeoffLengthx = PointToNumx(TakeoffLength)





    def designpoint(points):
        x_val = min(MinSpeedy[0], landingLengthy[0])
        i = 1
        y = searchXreturnY(points[0], x_val)
        pos = 0
        while (i<len(points)):
            curr = searchXreturnY(points[i], x_val)
            prev = searchXreturnY(points[pos], x_val)
            if(curr>prev):
                y = curr
                pos = i
            #print(y)
            i+=1
        point = Point(x_val, y*1.05)
        return point


    DesPoint = designpoint(Points)
    Sw = MassEstimate.Maxmass*9.81/DesPoint.x
    #print("THE WING AREA IS", Sw)
    b = math.sqrt(Sw*DragPolar.AR)
    Thrust = MassEstimate.Maxmass*9.81*DesPoint.y/1000
    #print("THE THRUST PER ENGINE IS", Thrust/2)


    class wingStuff:
        sweepQ = math.acos(1.16/(tlar.Mach+0.5))*(180/math.pi)
        taper = 0.2*(2-sweepQ*(math.pi/180))
        Sw = Sw
        b = b
        Thrust = Thrust
        Cr = 2*Sw/((1+taper)*b)
        Ct = Cr*taper
        MAClen = (2*Cr/3)*((1+taper+taper**2)/(1+taper))
        dih = 4
        


    class ClassII:
        kgs = 0.4535924
        meter = 0.3048
        lbs = 1/kgs
        ft = 1/meter
        gallon = 219.9692
        m_cube = 1/gallon

        # Global / common
        W_dg = MassEstimate.Maxmass*lbs   # design gross weight, lb
        N_z = 1.5*(2.1+(24000/(MassEstimate.Maxmass+10000)))  # ultimate load factor (1.5 * limit),
        q    = 0.5*rhoCruise*(vCruise**2)*(lbs/(ft**2))   # dynamic pressure at cruise, lb/ft^2

        # Wing (eq. 15.46)
        S_w      = Sw*(ft**2)  # wing planform area, ft^2
        W_fw     = 3167*lbs # fuel weight, lb
        A_w      = DragPolar.AR  # wing aspect ratio,
        Lambda_w = wingStuff.sweepQ  # wing sweep (typically @ 25% MAC), deg
        lambda_w = wingStuff.taper  # wing taper ratio,
        t_c_w    = 0.12  # wing thickness-to-chord ratio,

        # Horizontal tail (eq. 15.47)
        S_ht      = 10.33 *(ft**2) # horizontal tail area, ft^2
        A_ht      = 3.5  # horizontal tail aspect ratio,
        Lambda_ht = 25  # horizontal tail sweep (typically @ 25% MAC), deg
        lambda_ht = 0.7  # horizontal tail taper ratio,
        t_c_ht    = 0.12  # horizontal tail thickness-to-chord ratio,

        # Vertical tail (eq. 15.48)
        S_vt      = 10.25 *(ft**2) # vertical tail area, ft^2
        A_vt      = 1.7  # vertical tail aspect ratio,
        Lambda_vt = 35  # vertical tail sweep (typically @ 25% MAC), deg
        lambda_vt = 0.7  # vertical tail taper ratio,
        t_c_vt    = 0.12  # vertical tail thickness-to-chord ratio,
        H_t       = 4.17 *ft # horizontal tail height above fuselage, ft     
        H_v       = 4.17 *ft # vertical tail height above fuselage, ft         

        # Fuselage (eq. 15.49)
        S_f      = 70.27 *(ft**2)  # fuselage wetted area, ft^2
        L_t_fus  = 14.63 *ft # total fuselage length, ft
        L_over_D = 14.63/2.18  # fuselage fineness ratio (length/diameter),
        Vpr = 14.63*math.pi*1*(ft**3) #volume of pressurised area      #CHANGE THIS ONE
        Pdelta = 10.9 #cabin pressure
        W_press  = 11.9+(Vpr*Pdelta)**0.271 # pressurization weight allowance, lb

        # Landing gear (eqs. 15.5 - 15.51)
        N_l_main = 0.92*1.5  # landing load factor (Main),
        N_l_nose  = 0.15*1.5 # landing load factor (Nose), 
        W_l      = 0.91*MassEstimate.Maxmass*lbs # landing design gross weight, lb
        L_m_in   = 1.21 *ft*12 # main landing-gear strut length, in
        L_n_in   = 1.21 *ft*12 # nose landing-gear strut length, in

        # Installed engines (eq. 15.52)
        W_en  = 500 *lbs # single engine weight, lb
        N_en  = 2  # number of engines,
        # Fuel system (eq. 15.53)
        V_t   = 4.83 *gallon # total fuel volume, gal
        V_i   = V_t # integral/protected tanks volume (subset of V_t), gal
        N_t   = 2  # number of fuel tanks,

        # Flight controls (eq. 15.54)
        L_air = 14.63 *ft # aircraft total length, ft
        B_w   = 20.576 *ft # wingspan, ft

        # Hydraulics (eq. 15.55) - uses W_dg only

        # Electrical (eq. 15.56) uses W_fuel_system and W_avionics (outputs), listed here for completeness
        # (no extra inputs beyond those already defined)

        # Avionics (eq. 15.57)
        W_uav = 1100 # uninstalled avionics weight, lb

        # Air-conditioning & anti-ice (eq. 15.58)
        N_p = 10  # number of personnel onboard (crew + pax)
        M   = 0.68  # cruise Mach number

        # Furnishings (eq. 15.59) - uses W_dg only

        
    def compute_ga_class2_weights(p: ClassII) -> dict:

        def cosd(deg: float) -> float:
            return math.cos(math.radians(deg))

        # --- Wing (15.46)
        W_wing = (
            0.036
            * (p.S_w ** 0.758)
            * (p.W_fw ** 0.0035)
            * ((p.A_w / (cosd(p.Lambda_w) ** 2)) ** 0.6)
            * (p.q ** 0.006)
            * (p.lambda_w ** 0.04)
            * (((100.0 * p.t_c_w) / cosd(p.Lambda_w)) ** -0.3)
            * ((p.N_z * p.W_dg) ** 0.49)
        )

        # --- Horizontal tail (15.47)
        W_ht = (
            0.016
            * ((p.N_z * p.W_dg) ** 0.414)
            * (p.q ** 0.168)
            * (p.S_ht ** 0.896)
            * (((100.0 * p.t_c_ht) / cosd(p.Lambda_ht)) ** -0.12)
            * ((p.A_ht / (cosd(p.Lambda_ht) ** 2)) ** 0.043)
            * (p.lambda_ht ** -0.02)
        )

        # --- Vertical tail (15.48)
        Ht_over_Hv = (p.H_t / p.H_v) if getattr(p, "H_v", 0) else 0.0
        W_vt = (
            0.073
            * (1.0 + 0.2 * Ht_over_Hv)
            * ((p.N_z * p.W_dg) ** 0.376)
            * (p.q ** 0.122)
            * (p.S_vt ** 0.873)
            * (((100.0 * p.t_c_vt) / cosd(p.Lambda_vt)) ** -0.49)
            * ((p.A_vt / (cosd(p.Lambda_vt) ** 2)) ** 0.357)
            * (p.lambda_vt ** 0.039)
        )

        # --- Fuselage (15.49)
        W_fuselage = (
            0.052
            * (p.S_f ** 1.086)
            * ((p.N_z * p.W_dg) ** 0.177)
            * (p.L_t_fus ** -0.051)
            * (p.L_over_D ** -0.072)
            * (p.q ** 0.241)
            + p.W_press
        )

        # --- Landing gear (main & nose) (15.50–15.51)
        W_main_lg = 0.095 * ((p.N_l_main * p.W_l) ** 0.768) * ((p.L_m_in / 12.0) ** 0.409)
        W_nose_lg = 0.125 * ((p.N_l_nose * p.W_l) ** 0.566) * ((p.L_n_in / 12.0) ** 0.845)

        # --- Installed engines (total) (15.52)
        W_installed_engine = 2.575 * (p.W_en ** 0.922) * p.N_en

        # --- Fuel system (15.53)
        # Guard terms for Vi/Vt
        Vi_over_Vt = (p.V_i / p.V_t) if p.V_t else 0.0
        W_fuel_system = (
            2.49
            * (p.V_t ** 0.726)
            * ((1.0 / (1.0 + Vi_over_Vt)) ** 0.363)
            * (p.N_t ** 0.242)
            * (p.N_en ** 0.157)
        )

        # --- Flight controls (15.54)
        W_flight_controls = 0.053 * (p.L_air ** 1.536) * (p.B_w ** 0.371) * (
            (p.N_z * p.W_dg * 1.0e-4) ** 0.80
        )

        # --- Hydraulics (15.55)
        W_hydraulics = 0.001 * p.W_dg

        # --- Avionics (15.57)
        W_avionics = 2.117 * (p.W_uav ** 0.933)

        # --- Electrical (15.56)  (depends on W_fuel_system and W_avionics)
        W_electrical = 12.57 * ((W_fuel_system + W_avionics) ** 0.51)

        # --- Air conditioning & anti-ice (15.58)
        W_ac_antiice = 0.265 * (p.W_dg ** 0.52) * (p.N_p ** 0.68) * (W_avionics ** 0.17) * (p.M ** 0.08)

        # --- Furnishings (15.59)
        W_furnishings = 0.0582 * p.W_dg - 65.0

        return ((W_wing+W_ht+W_vt+W_fuselage+W_main_lg+W_nose_lg+W_installed_engine+W_fuel_system+W_flight_controls+W_hydraulics+W_electrical+W_avionics+W_ac_antiice+W_furnishings)*0.4535924+3167+1010)

        





    weights = compute_ga_class2_weights(ClassII)



    MassEstimate.Maxmass = weights

    #print(weights)





    #plt.plot(MinSpeedy, bruh, label='speed')
    #plt.plot(landingLengthy, bruh, label='Landing')
    #plt.plot(cruiseSpeedx, cruiseSpeedy, label='CruiseSpeed')
    #plt.plot(climbRatex, climbRatey, label='ClimbRate')
    #plt.plot(Cs119x, Cs119y, label='Cs119')
    #plt.plot(Cs121ax, Cs121ay, label='Cs121a')
    #plt.plot(Cs121bx, Cs121by, label='Cs121b')
    #plt.plot(Cs121cx, Cs121cy, label='Cs121c')
    #plt.plot(Cs121dx, Cs121dy, label='Cs121d')
    #plt.plot(TakeoffLengthx, TakeoffLengthy, label='Takeoff')


    #plt.plot(DesPoint.x, DesPoint.y,'mo', label='DesPoint')
    #plt.title('Matching Diagram')
    #plt.xlabel('Wing Loading')
    #plt.ylabel('Thrust to Weight')
    #plt.legend()

    #plt.show()
    print("THE WEIGHT IS ", weights)
    Cycle+=1

