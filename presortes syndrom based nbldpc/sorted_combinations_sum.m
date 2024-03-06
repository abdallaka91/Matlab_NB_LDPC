function result = sorted_combinations_sum(arr1, arr2)
    result = [];
    i = 1;
    j = 1;
    while i <= length(arr1) && j <= length(arr2)
        if arr1(i) <= arr2(j)
            result = [result, arr1(i) + arr2];
            i = i + 1;
        else
            result = [result, arr2(j) + arr1];
            j = j + 1;
        end
    end
    while i <= length(arr1)
        result = [result, arr1(i) + arr2];
        i = i + 1;
    end
    while j <= length(arr2)
        result = [result, arr2(j) + arr1];
        j = j + 1;
    end
    result = sort(result);
end


