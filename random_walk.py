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
def random_walk(steps: int) -> None:
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
    random_walk(STEPS) if i % 2 == 0 else random_walk(STEPS + 1)

x_locations = range(0, STEPS)
y_locations = range(0, STEPS)

heatmap_array = np.array(heatmap)


fig, ax = plt.subplots()
im = ax.imshow(heatmap_array, cmap = 'hot') # cmap = 'autumn'
ax.invert_yaxis()
plt.xticks(range(0, N))
plt.yticks(range(0, N))
plt.colorbar(im)
plt.title(f"Random Walk Simulation With {STEPS} Steps and {ITERATIONS} Iterations")
plt.show()

END_TIME = time.time()
print(f"Time Elapsed: {END_TIME - START_TIME}s")