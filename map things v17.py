import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point, MultiPoint
import random
import string
from shapely.ops import unary_union, cascaded_union

def generate_voronoi_diagram(width, height, num_cells, min_dist=5):
    points = []
    attempts = 0
    while len(points) < num_cells and attempts < 1000:
        point = np.random.rand(2) * [width, height]
        if all(np.linalg.norm(point - p) >= min_dist for p in points):
            points.append(point)
        attempts += 1
    return Voronoi(points)

def add_vertices_to_polygon(polygon, target_vertices, max_distance=10):
    # Ensure the target number of vertices is at least 30
    target_vertices = max(target_vertices, 30)

    # Simplify adding vertices to avoid overly complex shapes
    if len(polygon.exterior.coords) > target_vertices:
        return polygon.simplify(0.5, preserve_topology=False)
    else:
        new_coords = polygon.exterior.coords[:-1]  # Drop duplicate end point for processing
        while len(new_coords) < target_vertices:
            # Double the number of points by adding midpoints with random offsets
            new_coords = [new_coords[i // 2] if i % 2 == 0 else ((new_coords[i // 2][0] + new_coords[(i // 2 + 1) % len(new_coords)][0]) / 2 + random.uniform(-0.5, 0.5),
                                                                 (new_coords[i // 2][1] + new_coords[(i // 2 + 1) % len(new_coords)][1]) / 2 + random.uniform(-0.5, 0.5))
                          for i in range(2 * len(new_coords))]

            # Apply a maximum distance constraint
            new_coords = [new_coords[i] if np.linalg.norm(np.array(new_coords[i]) - np.array(new_coords[i - 1])) < max_distance
                          else ((new_coords[i][0] + new_coords[i - 1][0]) / 2, (new_coords[i][1] + new_coords[i - 1][1]) / 2)
                          for i in range(len(new_coords))]

        return Polygon(new_coords)

def generate_name():
    vowels = "aeiou"
    consonants = "".join(set(string.ascii_lowercase) - set(vowels))  # All consonants
    # Define less frequent consonants
    less_frequent_consonants = "zqx"
    more_frequent_consonants = "".join([c for c in consonants if c not in less_frequent_consonants])
    
    # Adjusting the length distribution
    if random.random() < 0.75:
        length = random.randint(5, 10)
    else:
        length = random.randint(3, 16)

    name_parts = []
    prev_char = ''
    prev_char_vowel = False
    for _ in range(length):
        if prev_char_vowel:  # If the previous character was a vowel
            # Decide on picking a consonant
            if random.random() < 0.8:  # Higher chance to pick a more frequent consonant
                new_char = random.choice(more_frequent_consonants)
            else:  # Significantly lower chance for 'z', 'x', 'q'
                new_char = random.choice(less_frequent_consonants)
        else:  # If the previous character was a consonant
            # Decide on picking a vowel or a more frequent consonant
            if random.random() < 0.6:  # Adjust probability to favor vowels after a consonant
                new_char = random.choice(vowels)
            else:
                new_char = random.choice(more_frequent_consonants)
        
        # Check for immediate repeated letters
        if new_char == prev_char:
            continue
        
        # Check for three or more of the same letter in a row
        if len(name_parts) >= 2 and new_char == name_parts[-1] == name_parts[-2]:
            continue

        name_parts.append(new_char)
        prev_char = new_char
        prev_char_vowel = new_char in vowels

    name = ''.join(name_parts).capitalize()

    if length > 4 and random.random() < 0.2:
        space_pos = random.randint(3, min(12, length - 1))
        name = name[:space_pos] + ' ' + name[space_pos:]

    return name

def plot_countries_with_capitals(countries, islands=None, background_color='skyblue'):
    if islands is None:
        islands = {}
    fig, ax = plt.subplots(figsize=(15, 15))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    plt.axis('off')

    country_colors = {}

    for i, country in enumerate(countries):
        color = np.random.rand(3,)
        # Ensure country color is not too close to the background color
        while np.linalg.norm(color - np.array([135, 206, 235]) / 255) < 0.1:
            color = np.random.rand(3,)
        country_colors[i] = color
        ax.fill(*country.exterior.xy, color=color, edgecolor='black')

        # Place capital randomly within the country
        capital_location = find_random_point_in_polygon(country)
        ax.plot(capital_location.x, capital_location.y, marker='*', color='black', markersize=9)  # Adjust marker size as needed

        # Display country name centered
        country_centroid = country.centroid
        country_name = generate_name()
        ax.text(country_centroid.x, country_centroid.y, country_name, ha='center', va='center',
                color='white', fontsize=9, fontweight='bold',
                bbox=dict(facecolor='black', edgecolor='none', pad=3.0))
    
    # Plot islands with the same color as their corresponding countries
    for country_index, island in islands.items():
        ax.fill(*island.exterior.xy, color=country_colors[country_index], edgecolor='black')

    plt.show()

def find_random_point_in_polygon(polygon):
    minx, miny, maxx, maxy = polygon.bounds
    while True:
        p = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
        if polygon.contains(p):
            return p

def create_and_plot_countries(num_countries=20, base_vertices=125):
    width, height = 100, 100
    vor = generate_voronoi_diagram(width, height, num_countries)
    countries = [Polygon(vor.vertices[region]) for region in vor.regions if -1 not in region and region != []]
    countries = [add_vertices_to_polygon(country, base_vertices * 3) for country in countries]

    # Randomly select 1-3 countries to generate islands for
    num_islands = random.randint(1, 3)
    islands = {}
    for i in range(num_islands):
        attempts = 0
        while attempts < 1000:
            # Generate an island with a random size and location
            island_width, island_height = random.uniform(5, 10), random.uniform(5, 10)
            island_location = np.random.rand(2) * np.array([width - island_width, height - island_height])
            island_points = [island_location + np.random.rand(2) * [island_width, island_height] for _ in range(random.randint(3, 6))]
            island = Polygon(island_points).convex_hull  # Use convex hull to ensure the polygon is valid

            # Check if island overlaps with any existing country
            overlap = False
            for country in countries:
                if island.intersects(country):
                    overlap = True
                    break

            if not overlap:
                islands[i] = island
                break
            attempts += 1

    plot_countries_with_capitals(countries, islands)

create_and_plot_countries()