import pandas as pd
import matplotlib.pyplot as plt
import argparse

def plot_data(df):
    # Normalize D and g values
    df['D_normalized'] = (df['D'] - df['D'].min()) / (df['D'].max() - df['D'].min())
    df['g_normalized'] = (df['g'] - df['g'].min()) / (df['g'].max() - df['g'].min())

    # Create the scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df['D_normalized'], df['g_normalized'], alpha=0.5)
    plt.xlabel('Normalized D')
    plt.ylabel('Normalized g')
    plt.title('Normalized D vs g')
    plt.grid(True)

    # Show the plot
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot D vs g from CSV file')
    parser.add_argument('csv_file', type=str, help='Path to the CSV file')
    args = parser.parse_args()

    # Read the CSV file
    df = pd.read_csv(args.csv_file)

    # Plot the data
    plot_data(df)

if __name__ == "__main__":
    main()
