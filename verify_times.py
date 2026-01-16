
import pandas as pd
from pathlib import Path
from nirsdexata.io import read_nirs
import datetime

def parse_excel_time(val):
    # val might be datetime.time object or string "MM:SS" or string "HH:MM"?
    # If pandas reads exact "18:00", it might be a string or time 18:00:00 (6 PM).
    # If the user means 18 minutes, it should likely be roughly 18.
    
    if isinstance(val, str):
        parts = val.split(':')
        if len(parts) == 2:
            return float(parts[0]) * 60 + float(parts[1])
    # If it is time object 18:00:00, that is 18 hours?? 
    # Or maybe it is 00:18:00?
    
    return val

def run():
    try:
        df = pd.read_excel("data-raw/participantdata.xlsx")
        # Ensure id is int
        df = df[pd.to_numeric(df['id'], errors='coerce').notnull()]
        df['id'] = df['id'].astype(int)
        
        print(f"{'ID':<5} {'Excel Time':<12} {'Duration':<10} {'Calc Start':<10} {'Diff':<10}")
        print("-" * 50)
        
        for _, row in df.iterrows():
            pid = row['id']
            time_val = row['time'] # Unparsed
            
            # Construct filename
            # Some are 'ptp1.txt', some 'PTP13_TSI.txt'.
            # Based on file list, mostly 'ptpX.txt'.
            fpath = Path(f"data-raw/ptp{pid}.txt")
            if not fpath.exists():
                # Check for other variants
                files = list(Path("data-raw").glob(f"ptp{pid}*.txt"))
                if not files:
                    files = list(Path("data-raw").glob(f"PTP{pid}*.txt"))
                
                if files:
                    fpath = files[0]
                else:
                    # print(f"File for {pid} not found")
                    continue
            
            # Read metadata (use our io module)
            try:
                # We only need header info, reading full file might be slow but safe
                rd = read_nirs(str(fpath))
                sfreq = rd.info.sfreq
                num_samples = len(rd.data)
                duration = num_samples / sfreq
                
                # Predict start
                # Logic: Start = Duration - 16*60 (960s)
                # If negative, 0
                calc_start = max(0, duration - 960)
                
                # Parse Excel time
                # Expecting format MM:SS?
                # If "18:00" -> 18 mins -> 1080s.
                # If we get 18:00:00 time object?
                excel_seconds = None
                if isinstance(time_val, datetime.time):
                     # Likely interpreted as HH:MM:SS. 
                     # If 18:00:00 -> 18 hours. Unlikely.
                     # If 00:18:00 -> 18 mins.
                     excel_seconds = time_val.hour * 3600 + time_val.minute * 60 + time_val.second
                elif isinstance(time_val, str):
                    parts = time_val.split(':')
                    if len(parts) >= 2:
                        excel_seconds = float(parts[0]) * 60 + float(parts[1])
                
                if excel_seconds is not None:
                    diff = abs(calc_start - excel_seconds)
                    print(f"{pid:<5} {str(time_val):<12} {duration:<10.2f} {calc_start:<10.2f} {diff:<10.2f}")
                else:
                    print(f"{pid:<5} {str(time_val):<12} {duration:<10.2f} {calc_start:<10.2f} ?")
                    
            except Exception as e:
                # print(f"Error reading {fpath}: {e}")
                pass

    except Exception as e:
        print(f"Main error: {e}")

if __name__ == "__main__":
    run()
