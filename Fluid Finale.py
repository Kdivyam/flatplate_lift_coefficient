#####################################################################
# File             : Fluid_Sim.py                                   #
# Function         : Computation of Lifting Flow about a Flat Plate #
# Author           : Divyam Khandelwal                              #
# Institute        : Mahindra Ecole Centrale                        #
#####################################################################

import math

import pylab as plt

x=[0 for i in range(59)]
y=[0 for j in range(60)]
h=[0 for i in range(59)]
k=[0 for j in range(60)]

#X-Coordinate Generation

for i in range(0,59):
    if(i<17):
        for i in range(0,17):
            x[i]=.5*(1-math.pow(1.1,(17-i)))
    elif(i<37):
        for i in range(17,37):
            x[i]=.05*(i-17)
    else:
        for i in range(37,59):
            x[i]=1.0-.5*(1-math.pow(1.1,(i-37)))

#Y-Coordinate Generation

for j in range(0,60):
    if(j<26):
        for j in range(0,26):
            y[j]=-.15+.5*(1-math.pow(1.1,(26-j)))
    elif(j<33):
        if(j<=29):
            for j in range(26,30):
                y[j]=.05*(j-29)
        elif(j==30):
            y[j]=y[j-1]
        else:
            y[j]=.05*(j-30)
    else:
        y[j]=.15-.5*(1-math.pow(1.1,(j-33)))

#H-Value Generation

for i in range(0,59):
    if(i>0 and i<58):
        h[i]=((x[i+1]-x[i-1])/2)

#K-Value Generation

for j in range(0,60):
    if(j==29):
        k[j]=(y[j+2]-y[j-1])/2
    elif((j>0 and j<29) or (j>30 and j<59)):
        k[j]=(y[j+1]-y[j-1])/2

#Initialising the solution matrix

phi_present=[[0 for j in range(59)] for i in range(60)]
phi_previous=[[0 for j in range(59)] for i in range(60)]

#Initialising the solution at points on grid

phi_upperplate=[0 for i in range(0,59)]
phi_lowerplate=[0 for i in range(0,59)]

#Specific Solution

v_inf=90
alpha=-4*math.pi/180
#l=[0 for i in range(58)]
#for i in range(37,58):
#    l[i]=(x[58]-x[i])/(x[58]-x[37])
l=[0 for i in range(58)]
for i in range(37,58):
    if(x[i]<1.06):
        l[i]=.75
    else:
        l[i]=0

w=1.7 #Weight
epsilon=.0001 #Error

state=False
counter=0

while(state==False):

    counter+=1

    for j in range (0,59):

        for i in range (0,58):

            phi_previous[j][i]=phi_present[j][i]

            if i==0 or i==58:
                pass

            elif j==0 or j==29 or j==30 or j==59:
                pass

            elif (i>16 and i<38) and j==31: #Upper, next-to plate NN11
                phi_present[j][i] =(((phi_present[j][i-1] + phi_present[j][i+1])*(k[j]**2)) + ((phi_present[j+1][i]+(k[j]*(v_inf*math.sin(alpha))))*(h[i]**2)))/((h[i]**2) + (2*(k[j]**2)))
                phi_present[j-1][i]=phi_present[j][i] +(k[j]*v_inf*math.sin(alpha))  #Updating upper plate point NN10

                phi_upperplate[i]=phi_present[j-1][i] #Storing required solution in separate list

            elif (i>16 and i<38) and j==28: #Lower, next-to plate NN13
                phi_present[j][i] =(((phi_present[j][i-1] + phi_present[j][i+1])*(k[j]**2)) + ((phi_present[j-1][i]-(k[j]*(v_inf*math.sin(alpha))))*(h[i]**2)))/((h[i]**2) + (2*(k[j]**2)))
                phi_present[j+1][i]=phi_present[j][i] -(k[j]*v_inf*math.sin(alpha)) #Updating upper plate point NN11

                phi_lowerplate[i]=phi_present[j+1][i] #Storing required solution in separate list

            elif (i>37 and i<58) and j==31: #Upper, next-to wake

                if i<57:
                    phi_present[j][i] =(((phi_present[j][i-1] + phi_present[j][i+1])*(k[j]**2)) + ((phi_present[j+1][i]+(k[j]*l[i]*(v_inf*math.sin(alpha))))*(h[i]**2)))/((h[i]**2) + (2*(k[j]**2))) #NN15
                    phi_present[j-1][i]=phi_present[j][i]+(k[j]*l[i]*v_inf*math.sin(alpha)) #Updating upper wake point NN19

                else:
                    phi_present[j][i] =(((phi_present[j][i-1])*(k[j]**2)) + ((phi_present[j+1][i]+(k[j]*l[i]*(v_inf*math.sin(alpha))))*(h[i]**2)))/((h[i]**2) + ((k[j]**2))) #NN17
                    phi_present[j-1][i]=phi_present[j][i] +(k[j]*l[i]*v_inf*math.sin(alpha)) #Updating upper wake point NN19

            elif (i>37 and i<58)and j==28: #Lower, next-to-wake
                if i<57:
                    phi_present[j][i] =(((phi_present[j][i-1] + phi_present[j][i+1])*(k[j]**2)) + ((phi_present[j-1][i]-(k[j]*l[i]*(v_inf*math.sin(alpha))))*(h[i]**2)))/((h[i]**2) + (2*(k[j]**2))) #NN16
                    phi_present[j+1][i]=phi_present[j][i] -(k[j]*l[i]*v_inf*math.sin(alpha)) #Updating lower wake point NN20

                else:
                    phi_present[j][i] =(((phi_present[j][i-1])*(k[j]**2)) + ((phi_present[j-1][i]-(k[j]*l[i]*(v_inf*math.sin(alpha))))*(h[i]**2)))/((h[i]**2) + ((k[j]**2))) #NN18
                    phi_present[j+1][i]=phi_present[j][i] -(k[j]*l[i]*v_inf*math.sin(alpha)) #Updating lower wake point NN20

            elif (i>0 and i<17) and j==29: #Lower side of backward extension

                if i>1:
                    phi_present[j][i] = (((phi_present[j][i+1]+phi_present[j][i-1])*(k[j]**2)) + ((phi_present[j+2][i]+phi_present[j-1][i])*(h[i]**2)))/(2*((h[i]**2) + (k[j]**2))) #NN21

                else:
                    phi_present[j][i] = (((phi_present[j][i+1])*(k[j]**2)) + ((phi_present[j+2][i]+phi_present[j-1][i])*(h[i]**2)))/((2*(h[i]**2)) + (k[j]**2)) #NN23

            elif (i>0 and i<17) and j==30: #Upper side of backward extension
                phi_present[j][i]=phi_present[j-1][i] #NN22


            elif i==1 and (j<28 or j>31):

                if j==1: #Bottom-left corner, next-to-boundary point NN2
                    phi_present[j][i] = ((phi_present[j][i+1]*(k[j]**2)) + (phi_present[j+1][i]*(h[i]**2)))/((h[i]**2) + (k[j]**2))

                elif j==58: #Top-left corner, next-to-boundary point NN3
                    phi_present[j][i] = ((phi_present[j][i+1]*(k[j]**2)) + (phi_present[j-1][i]*(h[i]**2)))/((h[i]**2) + (k[j]**2))

                else: #Next-to-left far field boundary line NN6
                    phi_present[j][i] = ((phi_present[j][i+1]*(k[j]**2)) + ((phi_present[j-1][i] + phi_present[j+1][i])*(h[i]**2)))/((2*(h[i]**2)) + (k[j]**2))

            elif i==57 and (j<28 or j>31):

                if j==1: #Bottom-right corner, next-to-boundary point NN5
                    phi_present[j][i]  = ((phi_present[j][i-1]*(k[j]**2)) + (phi_present[j+1][i]*h[i]**2))/((h[i]**2) + (k[j]**2))

                elif j==58: #Top-right corner, next-to-boundary point NN4
                    phi_present[j][i] = ((phi_present[j][i-1]*(k[j]**2)) + (phi_present[j-1][i]*(h[i]**2)))/((h[i]**2) + (k[j]**2))

                else: #Next-to-right far field boundary line NN8
                    phi_present[j][i] = ((phi_present[j][i-1]*(k[j]**2)) + ((phi_present[j+1][i] + phi_present[j-1][i])*(h[i]**2)))/((2*(h[i]**2)) + (k[j]**2))

            elif (i>1 and i<57) and j==1: #Next-to-lower far field boundary line NN9
                phi_present[j][i] = (((phi_present[j][i-1] + phi_present[j][i+1])*(k[j]**2)) + (phi_present[j+1][i]*(h[i]**2)))/((h[i]**2) + (2*(k[j]**2)))

            elif (i>1 and i<57) and j==58: #Next-to-upper far field boundary line NN7
                phi_present[j][i] = (((phi_present[j][i-1] + phi_present[j][i+1])*(k[j]**2)) + (phi_present[j-1][i]*(h[i]**2)))/((h[i]**2) + (2*(k[j]**2)))

            else: #Interior grid point NN1
                phi_present[j][i] = (((phi_present[j][i+1]+phi_present[j][i-1])*(k[j]**2)) + ((phi_present[j+1][i]+phi_present[j-1][i])*(h[i]**2)))/(2*((h[i]**2) + (k[j]**2)))

            phi_present[j][i] = w*phi_present[j][i]+((1-w)*phi_previous[j][i])

    #Updating solution at far-field boundary line points

    for i in (0,58):
        phi_present[0][i]=phi_present[1][i]
        phi_present[59][i]=phi_present[58][i]

    for j in (0,59):
        phi_present[j][0]=phi_present[j][1]
        phi_present[j][58]=phi_present[j][57]

    #Checking Convergence of Solution

    check=True

    for j in range(1,59):
        for i in range(1,58):
            if(abs(phi_previous[j][i]-phi_present[j][i])>epsilon):
                check=False
                break

    if(check==True):
        state=True

    print "Iteration:", counter

s = [[str(e) for e in row] for row in phi_present]
lens = [max(map(len, col)) for col in zip(*s)]
fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
table = [fmt.format(*row) for row in s]
print '\n'.join(table)


velocity_lowerplate=[0 for i in range(0,59)]
velocity_upperplate=[0 for i in range(0,59)]

for i in range(17,38):
    if(i==17):
        velocity_lowerplate[i]=(v_inf*math.cos(alpha))+((phi_lowerplate[i+1]-phi_lowerplate[i])/h[i])
    elif(i==37):
        velocity_lowerplate[i]=(v_inf*math.cos(alpha))+((phi_lowerplate[i]-phi_lowerplate[i-1])/h[i])
    else:
        velocity_lowerplate[i]=(v_inf*math.cos(alpha))+((phi_lowerplate[i+1]-phi_lowerplate[i-1])/(2*h[i]))


for i in range(17,38):
    if(i==17):
        velocity_upperplate[i]=(v_inf*math.cos(alpha))+((phi_upperplate[i+1]-phi_upperplate[i])/h[i])
    elif(i==37):
        velocity_upperplate[i]=(v_inf*math.cos(alpha))+((phi_upperplate[i]-phi_upperplate[i-1])/h[i])
    else:
        velocity_upperplate[i]=(v_inf*math.cos(alpha))+((phi_upperplate[i+1]-phi_upperplate[i-1])/(2*h[i]))


pressure_coeff_lower=[0 for i in range(0,59)]
pressure_coeff_upper=[0 for i in range(0,59)]
delta_c=[0 for i in range(0,59)]
pressure_coeff_lift=0

#Calculating pressure co-effecients

for i in range(17,38):
    pressure_coeff_upper[i]=1-math.pow((velocity_upperplate[i]/v_inf),2)
for i in range(17,38):
    pressure_coeff_lower[i]=1-math.pow((velocity_lowerplate[i]/v_inf),2)

# Calculating the lift

for i in range(17,38):
    pressure_coeff_lift+=(1/(1.05*(x[37]-x[17])))*((pressure_coeff_lower[i]-pressure_coeff_upper[i])*h[i])

for i in range(17,38):
    delta_c[i]=pressure_coeff_lower[i]-pressure_coeff_upper[i]

error=abs(pressure_coeff_lift-(2*math.pi*alpha))

print "Lift", pressure_coeff_lift
print "Error",error

#Plotting Velocity vs X-coordinate

vl_plot=velocity_lowerplate[17:38]
vu_plot=velocity_upperplate[17:38]
x_plot=x[17:38]
plt.plot(x_plot,vl_plot,linewidth=2.0,label='Lower plate GP')
plt.plot(x_plot,vu_plot,linewidth=2.0,label='Upper plate GP')
plt.legend(loc='upper right')
#pl.legend([plot1,plot2],('Upper plate grid points','Lower plate grid points'),'best')
plt.xlabel('X-Coordinate')
plt.ylabel('Velocity')
plt.title('Plot of Velocity vs. X ')
plt.show()

#Plotting Pressure Coeffecients vs X-coordinate

cu_plot=pressure_coeff_upper[17:38]
cl_plot=pressure_coeff_lower[17:38]
x_plot=x[17:38]
plt.plot(x_plot,cl_plot,linewidth=2.0,label='Lower plate CP')
plt.plot(x_plot,cu_plot,linewidth=2.0,label='Upper plate CP')
plt.legend(loc='upper right')
#pl.legend([plot1,plot2],('Upper plate grid points','Lower plate grid points'),'best')
plt.xlabel('X-Coordinate')
plt.ylabel('Pressure Coeffecient')
plt.title('Plot of Pressure Coeffecients vs. X ')
plt.show()

#Plotting Delta Cp vs X-coordinate

delta_c_plot=delta_c[17:38]
x_plot=x[17:38]
plt.plot(x_plot,delta_c_plot,linewidth=2.0)
plt.xlabel('X-Coordinate')
plt.ylabel(r'$\Delta$Cp')
plt.title('Plot of Delta Cp vs. X ')
plt.title(r'$\Delta$Cp vs. X')
plt.show()

"""
data = Data([
    Contour(
        z=[[10, 10.625, 12.5, 15.625, 20],
           [5.625, 6.25, 8.125, 11.25, 15.625],
           [2.5, 3.125, 5., 8.125, 12.5],
           [0.625, 1.25, 3.125, 6.25, 10.625],
           [0, 0.625, 2.5, 5.625, 10]]
    )
])
py.iplot(data)

for i in range(0,59):
    print i
    print phi_present[i][28]
for i in range(0,59):
    print i
    print phi_present[i][29]
for i in range(0,59):
    print i
    print phi_present[i][30]
for i in range(0,59):
    print i
    print phi_present[i][31]


cmap = mpl.colors.ListedColormap(['blue','black','red'])
bounds=[-6,-2,2,6]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# tell imshow about color map so that only set colors are used
img = pyplot.imshow(phi_present,interpolation='nearest',
                    cmap = cmap,norm=norm)

# make a color bar
pyplot.colorbar(img,cmap=cmap,
                norm=norm,boundaries=bounds,ticks=[-5,0,5])

pyplot.show()



    if((i>0 and i<16) and j==31): #Update equation for upper plate NN10
        phi_present[i][j-1]=phi_present[i][j]+(k[j]*v_inf*math.sin(alpha))

    if((i>0 and i<16)and j==28): #Update equation for lower plate NN12
        phi_present[i][j+1]=phi_present[i][j]-(k[j]*v_inf*math.sin(alpha))

    if((i>37 and i<58) and j==30): #Update equation for upper wake NN19
        phi_present[i][j-1]=(phi_present[i][j])+(k[j]*v_inf*l[i]*math.sin(alpha))

    if((i>37 and i<58) and j==29): #Update equation for lower wake NN20
        phi_present[i][j+1]=(phi_present[i][j])-(k[j]*v_inf*l[i]*math.sin(alpha))


#elif():#elif((i>0 and i<58) and (j>0 and j<59)): #Interior grid point NN1

            #    print "i: ",i
            #    print "j: ",j
            #    print "NN1"
            #    phi_present[i][j]=((math.pow(k[j],2)*(phi_present[i+1][j]+phi_present[i-1][j]))+(math.pow(h[i],2)*(phi_present[i][j+1]+phi_present[i][j-1])))/2*(math.pow(k[j],2)+math.pow(h[i],2))


    #phi_new[i][j]=(w*phi_present[i][j])+(1-w)*phi_present[i][j] #Play around with W?
    #phi_present[i][j]=phi_new[i][j]

    #Checking Convergence of Solution

data = Data([
    Contour(
        phi_present
    )
])
py.iplot(data)

#Updating equations

    for j in range(1,59):
        for i in range(1,58):
            if((i>16 and i<38) and j==30): #Update equation for upper plate NN10
                phi_present[i][j-1]=phi_present[i][j]+(k[j]*v_inf*math.sin(alpha))
                print "Needed1:",phi_present

            if((i>16 and i<38)and j==29): #Update equation for lower plate NN12
                phi_present[i][j+1]=phi_present[i][j]-(k[j]*v_inf*math.sin(alpha))
                print "Needed2:",phi_present

            if((i>37 and i<57) and j==30): #Update equation for upper wake NN19
                phi_present[i][j-1]=(phi_present[i][j])+(k[j]*v_inf*l[i]*math.sin(alpha))

            if((i>37 and i<57) and j==29): #Update equation for lower wake NN20
                phi_present[i][j+1]=(phi_present[i][j])-(k[j]*v_inf*l[i]*math.sin(alpha))

            #if((i>0 and i<17) and j==30): #Upper side of backward extension NN22
            #    phi_present[i][j]=(phi_present[i][j-1])
"""

