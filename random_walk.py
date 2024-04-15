# %%
import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time

# declare some constant variables
STEPS = 1000
ITERATIONS = 100000
N = 10 # coordinates will range from 0 - 99 x 0 - 99
START_TIME = time.time()

final_coordinates = []
heatmap = [[0]*N for _ in range(N)]  # create a heatmap of n x n size

# random walk function
def random_walk(steps):
    x = 0  # x init to 0
    y = 0  # y init to 0
    
    # iterate over "steps" amount of times
    for i in range(steps):
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
        # if x % n != x, reroll  
        # update coords
        if direction == "left":
            x -= 1
        elif direction == "right":
            x += 1
        elif direction == "up":
            y += 1
        elif direction == "down":
            y -= 1
        
    final_coordinates.append([x, y])
    heatmap[y][x] += 1

for i in range(ITERATIONS):
    random_walk(STEPS)

x_locations = range(0, STEPS)
y_locations = range(0, STEPS)

heatmap_array = np.array(heatmap)

fig, ax = plt.subplots()
im = ax.imshow(heatmap_array, cmap = 'hot')
plt.show()

END_TIME = time.time()
print(f"Time Elapsed: {END_TIME - START_TIME}s")

"""
# test visualizing plots

x_locations = range(0, STEPS)
y_locations = range(0, STEPS)

harvest = np.array([[0] * STEPS] * STEPS)

fig, ax = plt.subplots()
im = ax.imshow(harvest)

# Show all ticks and label them with the respective list entries
ax.set_xticks(np.arange(len(x_locations)), labels=x_locations)
ax.set_yticks(np.arange(len(y_locations)), labels=y_locations)

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.
for i in range(len(y_locations)):
    for j in range(len(x_locations)):
        text = ax.text(j, i, harvest[i, j],
                       ha="center", va="center", color="w")

ax.set_title("Harvest of local x_locations (in tons/year)")
#fig.tight_layout()
plt.show()"""