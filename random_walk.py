# %%
import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# declare some constant variables
STEPS = 10000
ITERATIONS = 10000

# random walk function
def random_walk(steps):
    x = 0  # x init to 0
    y = 0  # y init to 0
    
    #print(f"Initial position: ({x}, {y})")
    
    # iterate over "steps" amount of times
    for i in range(steps):
        direction = random.choice(["left", "right", "up", "down"])  # choose random direction
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
        
        #print(f"Step {i + 1}: ({x}, {y})")
    #print(f"Arrived at: ({x}, {y})")

for i in range(ITERATIONS):
    random_walk(STEPS)
print("done")
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