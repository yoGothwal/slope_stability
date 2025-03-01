import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scipy.integrate as spi

#parameters
cohesion = 10.00
friction_angle = 10

#slope geometry

crest_x = 10.000000000000000 
crest_y = 10.000000000000000
slope_angle = 45.00
slope_ratio = np.tan(np.radians(slope_angle))
slope_height =  crest_x * slope_ratio
  
# drawing boundaries of slope
x = [0, crest_x, crest_x + 50, crest_x + 50, -50, -50,0] 
y = [0, crest_x * slope_ratio, crest_x * slope_ratio, -50, -50,0,0]
plt.fill(x, y, color='blue', alpha=0.5, label="Slope")
plt.plot(x, y, 'ko-', linewidth="0.5") #slope's plot

# # drawing water
# wt_ratio = np.tan(np.radians(30))
# x = [0,crest_x +5 , crest_x + 50, crest_x + 50, -50, -50,0] 
# y = [0, (crest_x + 5)* wt_ratio, (crest_x) * ratio, -20, -20,0,0]
# plt.plot(x, y, '-bo')

# drawing failure circle
radius = 20 #using high precision
centre_x = 0 #circle x - cordinate
centre_y = 20 #circle y - cordinate

theta = np.linspace(0, 2 * np.pi,3000) # 3000 if good aprroximation.circle 
circle_x = radius * np.cos(theta) + centre_x
circle_y = radius * np.sin(theta) + centre_y # shifting the circle on y axis by 15 units

plt.plot(circle_x, circle_y, linestyle="--",color="k") #circle's plot
plt.fill(circle_x, circle_y, color='gray', alpha=0.5, label="Slope")

def equations(vars): #equations for finding cordinates of intersection
    x, y = vars
    eq1 = (x - centre_x) ** 2 + (y - centre_y) ** 2 - radius ** 2
    if x < 0:  
        eq2 = y  # Left boundary (horizontal ground, y = 0)
    elif x < crest_x:  
        eq2 = y - slope_ratio * x  # Sloped part of the surface
    else:  
        eq2 = y - crest_y  # Flat crest (top horizontal surface)

    return [eq1,eq2]# return true if the condition equals 0

ans1 = fsolve(equations,[-10,10])
ans2 = fsolve(equations, [20,crest_y])
x1, y1 = ans1 #intersection points from slope
x2, y2 = ans2 #intersection points from top slope side
x3, y3 = crest_x, crest_y  # Crest point
x4, y4 = 0, 0  # Toe point
print(f"Intersection Points:\n({x1:.2f}, {y1:.2f})\n({x2:.2f}, {y2:.2f})\nOther two points of slice:\n({x3:.2f}, {y3:.2f})\n({x4:.2f}, {y4:.2f})\n{0,0}")


def slice_area(x):
    circle_equation = centre_y - np.sqrt(radius ** 2 - (x - centre_x) ** 2)
    if x < 0:
        return 0 - circle_equation
    if x < crest_x:
        return slope_ratio * x - circle_equation
    else:
        return slope_height - circle_equation


total_area, _ = spi.quad(slice_area, x1,x2) 
print(f"total_area: {total_area:.4f}\n\n")

total_slices = 20
wedge_width = x2 - x1
slice_width = wedge_width / total_slices

def circle_y(x):
    return centre_y - np.sqrt(radius ** 2 - (x - centre_x) ** 2)

l = x2 - slice_width # left corner of slice
r = x2 #right corner
c = 0
sum = 0
while l >= x1:
    
    print("l value",l)
    c = c + 1
    circle_l = circle_y(l)
    circle_r = circle_y(r)
    mid = (l + r) / 2 #adding accuracy

    if l < crest_x:
        if r > crest_x:
            mid_y = (slope_height + (slope_ratio * l )) / 2
            y_coords = [circle_l, circle_r, slope_height, mid_y , slope_ratio * l]
        else:
            y_coords = [circle_l, circle_r, slope_ratio * r, slope_ratio * mid, slope_ratio * l]
        slice_height = slope_angle * l - circle_l
        
    else:
        slice_height = slope_height - circle_l
        y_coords = [circle_l, circle_r, slope_height, slope_height, slope_height] 
        
   
    x_coords = [l, r, r, mid, l]

    # Fill the slice
    plt.fill(x_coords, y_coords, color='red', edgecolor="black" )
    area, _ = spi.quad(slice_area, l, r)
    sum = sum + area
    print(f"Slice Number,{c} \nline: {l,r} \ncircle: {circle_l, circle_r} \n height:{slice_height} \narea: {area:.4f}\n")
    
    r = r - slice_width
    l = l - slice_width
    print(f"Slice {c}:")
    print(f"  x_coords: {x_coords}")
    print(f"  y_coords: {y_coords}")

print("Total sum after adding from each slice ",sum, "Total sum from integration: ", total_area)




plt.xlabel("X")
plt.ylabel("Y")
plt.title("Simple Circle")
plt.axis("equal")
plt.axhline(0, color='black', linewidth=1.5)
plt.axvline(0, color='black', linewidth=1.5)
plt.grid()
plt.show()