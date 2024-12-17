#figure 3E pie chart 

import matplotlib.pyplot as plt

# Sample data
labels = ['in-frame', 'out-of-frame']
sizes = [151, 122]  # Sizes for each slice, 151 is in frame and 122 out-of-frame
colors = ['black', 'dimgrey'] 

# Function to display actual values
def absolute_value(val):
    total = sum(sizes)
    absolute = int(val / 100. * total)
    return f'{absolute}'

# Create the pie chart
plt.figure(figsize=(8, 8))
plt.pie(sizes, labels=labels, colors=colors, autopct=absolute_value, startangle=140, textprops={'fontsize': 36})

# Add title
plt.title('FT with CDS/CDS')

# Display the plot
plt.show()
