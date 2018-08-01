from matplotlib import pyplot as plt
import sys

plt.pie((int(sys.argv[1]), int(sys.argv[2])), labels=("Enhancer", "No enhancer"))
plt.title("Relocalization peaks")
plt.savefig("relocalization_superenhancer_pie")
plt.close()
plt.pie((int(sys.argv[3]), int(sys.argv[4])), labels=("Enhancer", "No enhancer"))
plt.title("Background A compartment")
plt.savefig("background_superenhancer_pie")
