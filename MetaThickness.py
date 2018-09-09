import pylab as plt

coat_thickness = [0.2,1,2,4,8,16,32,48,64,80,96,112,128]

res_freq = [1.1462,1.1292,1.116,1.1036,1.0906,1.0838,1.08,1.0812,1.0847,1.0818,1.0714,1.0748,1.0782]

print "Thickness of Coating             Resonant Frequency"
for i in range(0,len(coat_thickness)):
    print "        ",coat_thickness[i],"                           ",res_freq[i]

plt.plot(coat_thickness , res_freq , linewidth=2.0)

plt.xlabel('Thickness of Coating ($\mu$m)')
plt.ylabel ('Resonant Frequency (THz)')

plt.title('Variation of Resonant Frequency with Thickness')
plt.show()