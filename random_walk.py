import random
import numpy as np
import matplotlib.pyplot as plt
import time

# declare some constant variables
STEPS = 10000
ITERATIONS = 1000000
N = 10 # coordinates will range from 0 - 9 x 0 - 9
START_TIME = time.time()

heatmap = [[0]*N for _ in range(N)]  # create a heatmap of n x n size

# random walk function
def random_walk(steps: int) -> None:
    x = 0  # x init to 0
    y = 0  # y init to 0
    # iterate over "steps" amount of times
    for _ in range(steps):
        # check possible moves at a location
        choices = []
        if x > 0:
            choices.append("left")
        if x < N - 1:
            choices.append("right")
        if y > 0:
            choices.append("down")
        if y < N - 1:
            choices.append("up")
        direction = random.choice(choices)  # choose random direction
        # update coords
        if direction == "left":
            x -= 1
        elif direction == "right":
            x += 1
        elif direction == "up":
            y += 1
        elif direction == "down":
            y -= 1
    # add iteration to heatmap
    heatmap[y][x] += 1

# call the random_walk function ITERATIONS number of times
for iteration in range(ITERATIONS):
    # alternate between even and odd steps
    random_walk(STEPS) if iteration % 2 == 0 else random_walk(STEPS + 1)

# convert heatmap 2d list to a numpy array for visualization purposes
heatmap_array = np.array(heatmap)

# visualize the heatmap with matplotlib
fig, ax = plt.subplots()
im = ax.imshow(heatmap_array, cmap = 'hot') # cmap = 'autumn'
# put origin on bottom left by inverting y ticks
ax.invert_yaxis()
# display all the numbers on the grid for the heatmap
plt.xticks(range(0, N))
plt.yticks(range(0, N))
# add a colorbar to the heatmap
plt.colorbar(im)
# add a title to the heatmap
plt.title(f"Random Walk Simulation With {STEPS} Steps and {ITERATIONS} Iterations")
# display heat map
plt.show()

for i in heatmap_array:
    for j in i:
        print(100*j/ITERATIONS)

# print the time elapsed
END_TIME = time.time()
print(f"Time Elapsed: {END_TIME - START_TIME}s")