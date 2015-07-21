cmap = lb.colormap()
cmap:append_color({1, 0, 0})
cmap:append_color({1, 1, 1})
cmap:append_color({0, 0, 1})

c = cmap:get_color(0)
print("color 0 = ",c[1],c[2],c[3])
c = cmap:get_color(1)
print("color 1 = ",c[1],c[2],c[3])
c = cmap:get_color(2)
print("color 2 = ",c[1],c[2],c[3])

print("Set color")

cmap:set_color(0,{0, 0, 1})
cmap:set_color(1,{1, 1, 1})
cmap:set_color(2,{1, 0, 0})

c = cmap:get_color(0)
print("color 0 = ",c[1],c[2],c[3])
c = cmap:get_color(1)
print("color 1 = ",c[1],c[2],c[3])
c = cmap:get_color(2)
print("color 2 = ",c[1],c[2],c[3])

print("Test assign value")
c = cmap:map_value(math.random())
print("Before\t",c[1],c[2],c[3])

c[1]=0
c[2]=0
c[3]=0

print("After\t", c[1],c[2],c[3])

print("Test write PNG")
nx=128
ny=128

rgb = lb.rgb(nx, ny)
for x = 0, nx - 1 do
	for y = 0, ny - 1 do
		c = cmap:map_value(math.random())
		rgb:set_pixel(x, ny - y - 1, c)
	end
end

c = rgb:get_pixel(nx/2,ny/2)
print("color = ",c[1],c[2],c[3])

rgb:save("filename.png")
collectgarbage() -- this should on C side

