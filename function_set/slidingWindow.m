function odindex = slidingWindow(data, windowSize, stepSize)
    em = cumtrapz(data.M.t, data.M.I.*data.M.U)/3600/1000; 
    em_d = [0 diff(em)];
    dataSize = size(data.M.t, 2); 
    numWindows = floor((dataSize - windowSize) / stepSize) + 1; 
    for i = 1:numWindows + 1
        startIdx = (i-1)*stepSize + 1;  
        endIdx = startIdx + windowSize - 1;  
        satart(i) = startIdx;
        endf(i) = endIdx;
    end
    endf(end) = min(endf(end), dataSize);
    em_sp = [];
    for i = 1:length(endf)
        em_sp = [em_sp em(endf(i))-em(satart(i))];
    end
    i = 1;
    n = length(satart);
    min_energy = 0.1;
    while i <= n
        if n == 1
            break
        end
        segment_energy = em_sp(i);
        if segment_energy < min_energy
            if i  == n
                satart(i) = [];
                endf(i - 1) = [];
                em_sp(i) = em_sp(i) + em_sp(i-1);
                em_sp(i) = [];
            else
                satart(i + 1) = [];
                endf(i) = [];
                em_sp(i) = em_sp(i) + em_sp(i+1);
                em_sp(i+1) = [];
            end
        else
            i = i + 1;
        end
        n = length(satart);
    end
    max_energy = 10;
    [satart, endf] = split_segments(em_d, satart, endf, max_energy);
    odindex = [satart; endf];
end

function [satart, endf] = split_segments(em_d, satart, endf, max_energy)
    i = 1;
    n = length(satart);
    while i <= n
        segment_energy = sum(em_d(satart(i):endf(i)));
        if segment_energy > max_energy
            split_points = find_split_points(em_d, satart(i), endf(i), max_energy);
            satart = [satart(1:i-1), split_points(1:end-1) + 1, satart(i+1:end)];
            endf = [endf(1:i-1), split_points(2:end), endf(i+1:end)];
            n = length(satart);
        else
            i = i + 1;
        end
    end
end

function split_points = find_split_points(em_d, start_idx, stop_idx, max_energy)
    n = 2;
    while true
        chunk_size = floor((stop_idx - start_idx) / n);
        split_indices = start_idx:chunk_size:stop_idx;
        energy = [];
        for ik = 1:length(split_indices) - 1
            energy = [energy, sum(em_d(split_indices(ik):split_indices(ik+1)))];
        end
        if energy < max_energy
            split_points = split_indices;
            break;
        else
            n = n+1;
        end
    end
    split_points(end) = stop_idx;
end