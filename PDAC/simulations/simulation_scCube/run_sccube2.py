import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.ndimage import gaussian_filter

# CONFIG
output_dir = "/home/martinpl/projects/datashare/grillus/art_visium_rect_custom"
os.makedirs(output_dir, exist_ok=True)

# SLIDE proportions: rectangular grid
n_rows = 60
n_cols = 60
spot_spacing = 1.0
spot_radius_mean = 0.48  # mean radius
spot_radius_jitter = 0.05  # amount of radius variation
position_jitter = 0.2  # jitter in position (to break perfect grid)

# Organic layers
n_layers = 4
n_colors = 4  # as per your image

# Custom colors: amarillo, verde, azul, gris
custom_colors = np.array([
    [252/255, 233/255, 79/255, 1.0],   # amarillo
    [138/255, 226/255, 52/255, 1.0],   # verde
    [52/255, 101/255, 164/255, 1.0],   # azul
    [186/255, 189/255, 182/255, 1.0]   # gris
])

# Loop for 50 images
for i in range(1, 51):
    print(f"Generating image {i}/50")

    xs = []
    ys = []

    # Build rectangular grid of spot centers + jitter
    for row in range(n_rows):
        for col in range(n_cols):
            x = col * spot_spacing + np.random.uniform(-position_jitter, position_jitter)
            y = row * spot_spacing + np.random.uniform(-position_jitter, position_jitter)
            xs.append(x)
            ys.append(y)

    xs = np.array(xs)
    ys = np.array(ys)

    # Create organic field
    grid_x = int(xs.max()) + 10
    grid_y = int(ys.max()) + 10
    field = np.zeros((grid_y, grid_x))

    # Superpose organic layers
    for _ in range(n_layers):
        cx = np.random.uniform(0, grid_x)
        cy = np.random.uniform(0, grid_y)
        amp = np.random.uniform(0.5, 1.5)
        sigma = np.random.uniform(10, 30)

        yy, xx = np.meshgrid(np.arange(grid_y), np.arange(grid_x), indexing='ij')
        layer = amp * np.exp(-((xx - cx)**2 + (yy - cy)**2) / (2 * sigma**2))
        field += layer

    # Smooth the field
    field = gaussian_filter(field, sigma=3)

    # Normalize
    field -= field.min()
    field /= field.max()

    # Sample field at spot positions
    spot_values = []
    for x, y in zip(xs, ys):
        xi = int(x)
        yi = int(y)
        val = field[yi, xi]
        spot_values.append(val)

    spot_values = np.array(spot_values)

    # Quantize to custom palette
    bins = np.linspace(0, 1, n_colors + 1)
    spot_bins = np.digitize(spot_values, bins) - 1

    # FIX → clip indices to avoid out-of-bounds
    spot_bins = np.clip(spot_bins, 0, n_colors - 1)

    # Map to custom colors
    colors_mapped = custom_colors[spot_bins]

    # PLOT
    fig, ax = plt.subplots(figsize=(8, 8))
    for x, y, color in zip(xs, ys, colors_mapped):
        # Vary radius slightly for organics
        radius = spot_radius_mean + np.random.uniform(-spot_radius_jitter, spot_radius_jitter)
        circle = plt.Circle((x, y), radius, color=color, ec='none')
        ax.add_patch(circle)

    ax.set_xlim(xs.min() - spot_spacing, xs.max() + spot_spacing)
    ax.set_ylim(ys.min() - spot_spacing, ys.max() + spot_spacing)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.tight_layout()

    # SAVE
    fig_path = os.path.join(output_dir, f"{i}.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()

print("✅ DONE →", output_dir)