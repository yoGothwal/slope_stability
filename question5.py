import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.integrate as spi
from matplotlib.patches import Arc

# Parameters
cohesion = 15.00
unit_weight = 18.00
saturated_unit_weight = 20.00
friction_angle = 10.00

# Drawing failure circle
radius = 15.00  # Using high precision
centre_x = 0.00  # Circle x-coordinate
centre_y = 15.00  # Circle y-coordinate

# Slope geometry
slope_angle = 45.00
wt_angle = 30.00
slope_ratio = np.tan(np.radians(slope_angle))
crest_x = 10.000000000000000
crest_y = slope_ratio * crest_x

slope_height = crest_x * slope_ratio
chord_mid = (1.00,1.00)
Rc = 5.00
reaction_angle = 45.00
X1 = X2 = 1.00
Y1 = Y2 = 1.00

# Function to calculate FOS_cohesion and FOS_friction for a given FOS
def calculate_fos(FOS):
    global chord_mid, Rc, reaction_angle, FOS_cohesion, FOS_friction, w_effective,Cm,P,Mw,Mc,Mp,arc_length,chord_length,alpha1, alpha, X1,X2,Y1,Y2

    # Drawing boundaries of slope
    x = [0, crest_x, crest_x + 50, crest_x + 50, -50, -50, 0]
    y = [0, crest_x * slope_ratio, crest_x * slope_ratio, -50, -50, 0, 0]
    plt.plot(x, y, 'ko-', linewidth="1", markersize=3)  # Slope's plot

    # Drawing water
    wt_ratio = np.tan(np.radians(30.00))
    x = [0, crest_x + 5, crest_x + 50, crest_x + 50, -50, -50, 0]
    y = [0, (crest_x + 5) * wt_ratio, (crest_x + 5) * wt_ratio, -20, -20, 0, 0]
    plt.plot(x, y, '-bo', markersize=3,linewidth="0.5")

    

    theta = np.linspace(0, 2 * np.pi, 3000)
    circle_x = radius * np.cos(theta) + centre_x
    circle_y = radius * np.sin(theta) + centre_y
    plt.plot(circle_x, circle_y, linestyle="--", color="red", linewidth="0.8")  # failure surface. Circle's plot

    def slope_boundary(x, y):
        if x < 0:
            return y  # Horizontal ground, y = 0
        elif x < crest_x:
            return y - slope_ratio * x  # Slope
        else:
            return y - crest_y  # Flat crest

    def circle_boundary(x, y):
        return (x - centre_x) ** 2 + (y - centre_y) ** 2 - radius ** 2

    def equations(vars):  # Equations for finding coordinates of intersection
        x, y = vars
        eq1 = circle_boundary(x, y)
        eq2 = slope_boundary(x, y)
        return [eq1, eq2]

    ans1 = fsolve(equations, [-10, 10])
    ans2 = fsolve(equations, [20, crest_y])
    X1, Y1 = ans1  # Left intersection points
    X2, Y2 = ans2  # Right intersection points

    def inclination_angle(coord1, coord2, ref_point):
        x1, y1 = coord1
        x2, y2 = coord2
        h, k = ref_point

        angle_rad = np.arctan2(y2 - k, x2 - h) - np.arctan2(y1 - k, x1 - h)
        angle_deg = np.degrees(angle_rad)

        return abs(angle_deg)

    coord1 = (X1, Y1)
    coord2 = (X2, Y2)
    ref_point = (centre_x, centre_y)
    alpha = inclination_angle(coord1, coord2, ref_point)

    def slope_equation(x):
        circle_equation = centre_y - np.sqrt(radius ** 2 - (x - centre_x) ** 2)
        if x < 0:
            return 0 - circle_equation
        if x < crest_x:
            return slope_ratio * x - circle_equation
        else:
            return slope_height - circle_equation

    def water_table_equation(x):
        circle_equation = centre_y - np.sqrt(radius ** 2 - (x - centre_x) ** 2)
        if x < 0:
            return 0 - circle_equation
        if x < crest_x + 5:
            return wt_ratio * x - circle_equation
        else:
            return wt_ratio * crest_x - circle_equation

    def slice_area(x):  # Gives equation of circular failure surface
        return slope_equation(x)

    def slice_area_wet(x):  # Gives equation of circular failure surface
        return water_table_equation(x)

    # Total area from x1 to x2
    area, _ = spi.quad(slice_area, X1, X2)
    # Wet area from x1 to x2
    area_wet, _ = spi.quad(slice_area_wet, X1, X2)
    area_dry = area - area_wet

    # Calculating effective weight
    w_effective = area_dry * unit_weight + area_wet * np.cos(np.radians(wt_angle)) * (saturated_unit_weight - unit_weight)

    # Driving moment
    coord1 = (centre_x, 0) # coord1 is choosen on y axis of circle. 
    chord_mid = ((X1 + X2) / 2, (Y1 + Y2) / 2)  # Chord middle point
    ref_point = (centre_x, centre_y) #point where angle is subtended
    alpha1 = inclination_angle(coord1, chord_mid, ref_point)  # Angle which W makes with y axis between coord1, ref_point, coord2


    lever = radius * np.sin(np.radians(alpha1))
    Mw = lever * w_effective
    # print("alpha1", alpha1)

    # Cohesion
    def chord_length(angle):
        return 2 * radius * np.sin(np.radians(alpha) / 2)

    def arc_length(angle):
        return (2 * np.pi * radius) * alpha / 360

    Rc = (arc_length(alpha) / chord_length(alpha)) * radius  # Slightly more than circle radius
    reaction_angle = alpha1 - friction_angle
    # print("Arc length:\n", arc_length(alpha))
    # print("Chord length:\n", chord_length(alpha))
    # print("Radius of chord, Rc\n", Rc)
    # print("Angle between W and Cm\n", 90 + alpha1)
    # print("Angle between reaction P and vertical axis:\n", reaction_angle)

    Cm = (cohesion / FOS) * arc_length(alpha)  # Total cohesive force along the arc
    cohesion_mobilized = Cm / chord_length(alpha)
    FOS_cohesion = cohesion / cohesion_mobilized

    Mc = Rc * Cm

    # Reaction
    P = np.sqrt(w_effective * w_effective + Cm * Cm + 2 * w_effective * Cm * np.cos(np.radians(90.00 + alpha1)))  # Resultant of Cm and W

    Mp = P * (radius * np.sin(np.radians(friction_angle)))
    # FOS = 1 => Driving moment = Resisting moment
    FOS_friction = (Mp + Mc) / Mw
    
    plt.axhline(0, color='black', linewidth=0.8 )
    plt.axvline(0, color='black', linewidth=0.8 )
    plt.xlabel("x(m)")
    plt.ylabel("y(m)")
    print("FOS assumed: ", FOS)
    print("New FOS (cohesion): ", FOS_cohesion)
    print("New FOS (phi): ", FOS_friction)
    


    return FOS_cohesion, FOS_friction

# Generating  FOS values
FOS_values = np.concatenate((np.linspace(0.01,1,101),np.linspace(1,10,1001)))  # Example FOS values
print(FOS_values)
FOS_assumed = 10
minimum = 10000
FOS_cohesion_values = []
FOS_friction_values = []

for FOS in FOS_values:
    FOS_cohesion, FOS_friction = calculate_fos(FOS)
    diff = np.abs(FOS_cohesion - FOS_friction)
    print("del FOS", diff)
    print("\n")

    if diff < minimum:
        minimum = diff
        FOS_assumed = FOS
    
    FOS_cohesion_values.append(FOS_cohesion)
    FOS_friction_values.append(FOS_friction)

def plot_point(x, y):
    plt.scatter(x, y, color='black', marker='o', s=10, label='Center')
    plt.text(x, y, f'({x:.2f}, {y:.2f})', fontsize=8, verticalalignment='bottom', horizontalalignment='right')
plot_point(centre_x, centre_y)
# Drawing boundaries of chord
plot_point(X1, Y1)
plot_point(X2, Y2)
plt.plot((X1, X2), (Y1, Y2), color="k", linewidth="0.5")

##########################
# calculating coordinates of location of Rc
def line(x, y, coord):
    x_val, y_val = coord
    return y - ((y_val - centre_y) / (x_val - centre_x)) * (x - x_val) - y_val

def cohesion_boundary(x, y):
    return (x - centre_x) ** 2 + (y - centre_y) ** 2 - Rc ** 2

def equations(vars):  # Equations for finding coordinates of intersection
    x, y = vars
    eq1 = line(x, y, chord_mid )
    eq2 = cohesion_boundary(x, y)
    return [eq1, eq2]

ans1 = fsolve(equations, (9.00,2.00))
ans2 = fsolve(equations, (-8.00,29.00))
x1, y1 = ans1 # coordinates where Rc is active
x2, y2 = ans2
x3, y3 = chord_mid
plt.scatter(x1, y1, color='red', s=5)  # drawing coordinates where Rc is active
print("\nPoint of application of Rc or P",(x1, y1), (x2, y2))

##################################

print(f"Intersection Points:\n({X1:.2f}, {Y1:.2f})\n({X2:.2f}, {Y2:.2f})\n")
print("Arc length:\n", arc_length(alpha))
print("Chord length:\n", chord_length(alpha))
print("Radius of chord, Rc\n", Rc)
print(f"\nW:{w_effective},\n Cm:{Cm},\n P:{P},\n Mw:{Mw},\n Mc:{Mc},\n Mp:{Mp},\n")
print("Angle between W and Cm\n", 90.00 + alpha1)
print("Angle between reaction P and vertical axis :\n", reaction_angle)


def draw_circle(x, y, r):
    theta = np.linspace(0, 2 * np.pi, 300)
    plt.plot(x + r * np.cos(theta), y + r * np.sin(theta), 'k-', linewidth=0.5)

def draw_line(coord1, coord2):
    x , y = zip(coord1, coord2)
    plt.plot(x, y, "k--", linewidth="0.5") 

def draw_line_with_slope(coord1, m):
    x1, y1 = coord1
    x = np.linspace(x1 -  radius, x1 +  radius/2,11)
    y = m * (x - x1) + y1
    plt.plot(x, y, "k--",linewidth="0.5")

draw_line((x1, y1), (centre_x, centre_y))
draw_line_with_slope((x1,y1),  np.tan(np.radians(90.00 + reaction_angle))) 
draw_circle(centre_x, centre_y , Rc)#circle where Cm is
draw_circle(centre_x, centre_y , Rc * np.sin(np.radians(friction_angle)) )# drawing friction circle
print("FOS assumed: ", FOS_assumed)
print("New FOS (cohesion): ", FOS_cohesion)
print("New FOS (phi): ", FOS_friction)

plt.figure(figsize=(8, 8))  # Ensuring a square figure to match aspect ratio
plt.plot(FOS_cohesion_values, FOS_friction_values, 'bo-', label='FOS_cohesion vs FOS_friction',linewidth="1")
draw_line((-1.00,-1.00), (10.00,10.00))


plt.axhline(0, color='black', linewidth=0.8, )
plt.axvline(0, color='black', linewidth=0.8, ) 

plt.xlabel("FOS_cohesion")
plt.ylabel("FOS_friction")
plt.title("FOS_cohesion vs FOS_friction")

plt.grid()
plt.legend()

# Ensure the aspect ratio is equal
plt.gca().set_aspect('equal', adjustable='datalim')

plt.show()