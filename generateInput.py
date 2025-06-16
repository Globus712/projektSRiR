import random
import sys

size = 100
adj_list = []
propability = 0.85
max_weight = 200

if __name__ == "__main__":
    if len(sys.argv) > 1:
        size = int(sys.argv[1])

    if len(sys.argv) > 2:
        propability = float(sys.argv[2])

    for i in range(size):
        adj_list.append([])

        for j in range(size):
            random_number = random.random()

            if random_number < 0.9:
                random_weight = random.randint(5, 200)
                adj_list[i].append((j, random_weight))

    with open(f"./input/input_{size}.txt", "w") as file:
        file.write(f"{size}\n")
        
        for x in adj_list:
            file.write(f"{len(x)}")

            for y in x:
                file.write(f" {y[0]},{y[1]}")
            
            file.write("\n")
    # for x in adj_list:
    #     for y in x:
    #         print(f"{y}, ", end=" ")
        
    #     print()