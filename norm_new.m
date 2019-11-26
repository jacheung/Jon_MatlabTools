function normalized = norm_new(array)
        
    max_response = max(array);
    normalized = array./ max_response;