
# %%
import random

# random walk function
def random_walk(steps):
    x = 0  # x init to 0
    y = 0  # y init to 0
    
    print(f"Initial position: ({x}, {y})")
    
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
        
        print(f"Step {i + 1}: ({x}, {y})")
        
        # stop running if arrive back to start
        if x == 0 and y == 0:
            print(f"Arrived back to start: ({x}, {y})")
            break

def main():
    while input:
        steps = int(input("Enter the number of steps for the random walk: "))
        random_walk(steps)

main()
