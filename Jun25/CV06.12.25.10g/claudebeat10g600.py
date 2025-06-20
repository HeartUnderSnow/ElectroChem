# -*- coding: utf-8 -*-
"""
Enhanced Cyclic Voltammetry Analysis for FLiNaK System
Includes all CSV files and improved analysis capabilities

@author: myheartUnderSnow
"""

import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np
from scipy.signal import find_peaks
import os
from pathlib import Path

# Define all file paths - combining original and new files
file_info = [
    # Original files from your script
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1617_CV_FLiNaK_Ni_601_200mVperS.csv', 'temp': 601, 'scan_rate': 200, 'time': '1617'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1623_CV_FLiNaK_Ni_601_200mVperS.csv', 'temp': 601, 'scan_rate': 200, 'time': '1623'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1659_CV_FLiNaK_Ni_601_200mVperS.csv', 'temp': 601, 'scan_rate': 200, 'time': '1659'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1701_CV_FLiNaK_Ni_600_200mVperS.csv', 'temp': 600, 'scan_rate': 200, 'time': '1701'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1703_CV_FLiNaK_Ni_600_200mVperS.csv', 'temp': 600, 'scan_rate': 200, 'time': '1703'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1723_CV_FLiNaK_Ni_600_100mVperS.csv', 'temp': 600, 'scan_rate': 100, 'time': '1723'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1723_CV_FLiNaK_Ni_600_200mVperS.csv', 'temp': 600, 'scan_rate': 200, 'time': '1723b'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1730_CV_FLiNaK_Ni_600_100mVperS.csv', 'temp': 600, 'scan_rate': 100, 'time': '1730'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1733_CV_FLiNaK_Ni_600_100mVperS.csv', 'temp': 600, 'scan_rate': 100, 'time': '1733'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1735_CV_FLiNaK_Ni_600_100mVperS.csv', 'temp': 600, 'scan_rate': 100, 'time': '1735'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1737_CV_FLiNaK_Ni_600_50mVperS.csv', 'temp': 600, 'scan_rate': 50, 'time': '1737'},
    {'path': '/home/galriel/Public/ElectroChem/Jun25/CV06.12.25.10g/1741_CV_FLiNaK_Ni_600_50mVperS.csv', 'temp': 600, 'scan_rate': 50, 'time': '1741'},
    
    # New files from uploaded CSVs
    {'path': '1453_CV_FLiNaK_Ni_600_100mVperS.csv', 'temp': 600, 'scan_rate': 100, 'time': '1453'},
    {'path': '1501_CV_FLiNaK_Ni_602_100mVperS.csv', 'temp': 602, 'scan_rate': 100, 'time': '1501'},
    {'path': '1503_CV_FLiNaK_Ni_602_100mVperS.csv', 'temp': 602, 'scan_rate': 100, 'time': '1503'},
    {'path': '1506_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1506'},
    {'path': '1508_CV_FLiNaK_Ni_602_100mVperS.csv', 'temp': 602, 'scan_rate': 100, 'time': '1508'},
    {'path': '1512_CV_FLiNaK_Ni_602_100mVperS.csv', 'temp': 602, 'scan_rate': 100, 'time': '1512'},
    {'path': '1516_CV_FLiNaK_Ni_602_50mVperS.csv', 'temp': 602, 'scan_rate': 50, 'time': '1516'},
    {'path': '1520_CV_FLiNaK_Ni_603_50mVperS.csv', 'temp': 603, 'scan_rate': 50, 'time': '1520'},
    {'path': '1523_CV_FLiNaK_Ni_603_50mVperS.csv', 'temp': 603, 'scan_rate': 50, 'time': '1523'},
    {'path': '1530_CV_FLiNaK_Ni_603_50mVperS.csv', 'temp': 603, 'scan_rate': 50, 'time': '1530'},
    {'path': '1537_CV_FLiNaK_Ni_602_50mVperS.csv', 'temp': 602, 'scan_rate': 50, 'time': '1537'},
    {'path': '1543_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1543'},
    {'path': '1545_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1545'},
    {'path': '1547_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1547'},
    {'path': '1548_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1548'},
    {'path': '1551_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1551'},
    {'path': '1553_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1553'},
    {'path': '1610_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1610'},
    {'path': '1612_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1612'},
    {'path': '1614_CV_FLiNaK_Ni_602_200mVperS.csv', 'temp': 602, 'scan_rate': 200, 'time': '1614'},
]

def load_cv_data(file_info_list):
    """Load all CV data files into a dictionary"""
    cv_data = {}
    
    for info in file_info_list:
        try:
            if os.path.exists(info['path']):
                df = pd.read_csv(info['path'])
                cv_data[info['time']] = {
                    'data': df,
                    'voltage': df['Potential (V)'],
                    'current': df['Current (A)'],
                    'temp': info['temp'],
                    'scan_rate': info['scan_rate'],
                    'time': info['time']
                }
                print(f"Loaded: {info['time']} - {info['temp']}°C - {info['scan_rate']}mV/s")
            else:
                print(f"File not found: {info['path']}")
        except Exception as e:
            print(f"Error loading {info['path']}: {e}")
    
    return cv_data

def find_k_redox_peaks(voltage, current, v_min, v_max):
    """Find K+/K redox peaks - anodic maxima and cathodic onset"""
    # Create mask for the potential range
    mask = (voltage >= v_min) & (voltage <= v_max)
    v_range = voltage[mask]
    i_range = current[mask]
    
    if len(v_range) == 0:
        return None, None
    
    original_indices = np.where(mask)[0]
    
    # Convert to numpy arrays for easier indexing
    i_array = i_range.values
    v_array = v_range.values
    
    # ANODIC PEAK: Find the maximum positive current peak
    positive_threshold = 0.002  # 2 mA threshold for anodic peaks
    positive_mask = i_array > positive_threshold
    
    anodic_index = None
    if np.any(positive_mask):
        pos_currents = i_array[positive_mask]
        pos_original_indices = original_indices[positive_mask]
        
        # Find anodic peak (highest positive current)
        max_current_idx = np.argmax(pos_currents)
        anodic_index = pos_original_indices[max_current_idx]
    
    # CATHODIC PEAK: Find the point closest to zero but still negative
    # in the very specific cathodic region
    cathodic_voltage_min = -1.69
    cathodic_voltage_max = -1.675
    
    # Filter for negative current in this narrow cathodic region
    cathodic_mask = (v_array >= cathodic_voltage_min) & (v_array <= cathodic_voltage_max) & \
                   (i_array < 0) & (i_array > -0.006)  # Small negative current only
    
    cathodic_index = None
    if np.any(cathodic_mask):
        cathodic_currents = i_array[cathodic_mask]
        cathodic_voltages = v_array[cathodic_mask]
        cathodic_original_indices = original_indices[cathodic_mask]
        
        # Find the point closest to zero current (least negative)
        closest_to_zero_idx = np.argmax(cathodic_currents)  # Maximum of negative values = closest to zero
        cathodic_index = cathodic_original_indices[closest_to_zero_idx]
    
    return anodic_index, cathodic_index

def calculate_redox_potential(voltage, anodic_idx, cathodic_idx):
    """Calculate redox potential from anodic and cathodic peak indices"""
    if anodic_idx is not None and cathodic_idx is not None:
        V_anodic = voltage.iloc[anodic_idx]
        V_cathodic = voltage.iloc[cathodic_idx]
        return (V_anodic + V_cathodic) / 2
    else:
        return None

def analyze_all_data(cv_data):
    """Analyze K+/K redox peaks for all datasets"""
    # Define potential range for K+/K redox
    k_redox_min = -1.75
    k_redox_max = -1.5
    
    results = {}
    
    print("K+/K Redox Potential Analysis:")
    print("="*70)
    
    for time_key, data in cv_data.items():
        voltage = data['voltage']
        current = data['current']
        temp = data['temp']
        scan_rate = data['scan_rate']
        
        # Find K+/K redox peaks
        anodic_idx, cathodic_idx = find_k_redox_peaks(voltage, current, k_redox_min, k_redox_max)
        
        # Calculate potentials and currents
        V_anodic = voltage.iloc[anodic_idx] if anodic_idx is not None else None
        V_cathodic = voltage.iloc[cathodic_idx] if cathodic_idx is not None else None
        I_anodic = current.iloc[anodic_idx] if anodic_idx is not None else None
        I_cathodic = current.iloc[cathodic_idx] if cathodic_idx is not None else None
        
        # Calculate redox potential
        E_redox = calculate_redox_potential(voltage, anodic_idx, cathodic_idx)
        
        # Store results
        results[time_key] = {
            'temp': temp,
            'scan_rate': scan_rate,
            'anodic_potential': V_anodic,
            'cathodic_potential': V_cathodic,
            'anodic_current': I_anodic,
            'cathodic_current': I_cathodic,
            'redox_potential': E_redox,
            'peak_separation': abs(V_anodic - V_cathodic)*1000 if (V_anodic is not None and V_cathodic is not None) else None
        }
        
        # Print results
        print(f"Time {time_key} - {temp}°C - {scan_rate}mV/s:")
        print(f"  Anodic Peak:      {V_anodic:.3f} V ({I_anodic:.6f} A)" if V_anodic is not None else "  Anodic Peak:      Not found")
        print(f"  Cathodic Peak:    {V_cathodic:.3f} V ({I_cathodic:.6f} A)" if V_cathodic is not None else "  Cathodic Peak:    Not found")
        print(f"  Redox Potential:  {E_redox:.3f} V" if E_redox is not None else "  Redox Potential:  Cannot calculate")
        if results[time_key]['peak_separation'] is not None:
            print(f"  Peak Separation:  {results[time_key]['peak_separation']:.1f} mV")
        print()
    
    return results

def plot_by_scan_rate(cv_data, scan_rates_to_plot=[50, 100, 200]):
    """Create separate plots for each scan rate"""
    
    for scan_rate in scan_rates_to_plot:
        plt.figure(figsize=(8, 6))
        
        # Filter data for this scan rate
        filtered_data = {k: v for k, v in cv_data.items() if v['scan_rate'] == scan_rate}
        
        if not filtered_data:
            continue
            
        # Plot each dataset for this scan rate
        for time_key, data in filtered_data.items():
            label = f"{time_key} ({data['temp']}°C)"
            plt.plot(data['voltage'], data['current'], label=label, alpha=0.8)
        
        plt.xlabel('Potential (V) vs. Ni', fontsize=14)
        plt.ylabel('Current (A)', fontsize=14)
        plt.title(f'CVs of FLiNaK at {scan_rate} mV/s', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show()

def plot_scan_rate_comparison(cv_data, temp_filter=None):
    """Plot comparison of different scan rates for a specific temperature"""
    plt.figure(figsize=(10, 7))
    
    # Group data by scan rate
    scan_rate_groups = {}
    for time_key, data in cv_data.items():
        if temp_filter is None or data['temp'] == temp_filter:
            scan_rate = data['scan_rate']
            if scan_rate not in scan_rate_groups:
                scan_rate_groups[scan_rate] = []
            scan_rate_groups[scan_rate].append((time_key, data))
    
    # Plot representative curves for each scan rate
    colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
    color_idx = 0
    
    for scan_rate in sorted(scan_rate_groups.keys()):
        # Take the first dataset for each scan rate as representative
        time_key, data = scan_rate_groups[scan_rate][0]
        color = colors[color_idx % len(colors)]
        
        plt.plot(data['voltage'], data['current'], 
                label=f'{scan_rate} mV/s ({data["temp"]}°C)', 
                color=color, linewidth=2)
        color_idx += 1
    
    plt.xlabel('Potential (V) vs. Ni', fontsize=14)
    plt.ylabel('Current (A)', fontsize=14)
    temp_str = f'at {temp_filter}°C' if temp_filter else 'All Temperatures'
    plt.title(f'CV Scan Rate Comparison {temp_str}', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

# Main execution
if __name__ == "__main__":
    # Load all CV data
    cv_data = load_cv_data(file_info)
    
    # Analyze all data for K+/K redox peaks
    results = analyze_all_data(cv_data)
    
    # Create summary statistics
    valid_redox_potentials = [r['redox_potential'] for r in results.values() if r['redox_potential'] is not None]
    if valid_redox_potentials:
        print(f"\nSummary Statistics for K+/K Redox Potentials:")
        print(f"Average: {np.mean(valid_redox_potentials):.3f} V")
        print(f"Std Dev: {np.std(valid_redox_potentials):.3f} V")
        print(f"Range: {np.min(valid_redox_potentials):.3f} to {np.max(valid_redox_potentials):.3f} V")
    
    # Create plots
    print("\nGenerating plots...")
    
    # Plot by scan rate
    plot_by_scan_rate(cv_data, [50, 100, 200])
    
    # Plot scan rate comparison for 600°C (if available)
    plot_scan_rate_comparison(cv_data, temp_filter=600)
    
    # Plot scan rate comparison for 602°C (if available)
    plot_scan_rate_comparison(cv_data, temp_filter=602)
    
    print("Analysis complete!")
