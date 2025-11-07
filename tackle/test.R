infer_ordered_factor <- function(column) { # this probably doesn't do what you want it to do
  # Function to infer ordering for a given column vector

  # Extract numeric values from strings
  numeric_values <- as.numeric(gsub("[^0-9.-]+", "", column))

  # Find the minimum numeric value
  min_num <- min(numeric_values, na.rm = TRUE)
  if (is.infinite(min_num)) {
    min_num <- 0  # Default to 0 if no numeric values are found
  }

  # Initialize order values with numeric values
  order_values <- numeric_values

  # Indices of non-numeric values
  non_numeric_indices <- which(is.na(order_values))

  if (length(non_numeric_indices) > 0) {
    non_numeric_values <- tolower(column[non_numeric_indices])

    # Assign special order value for "veh" or "vehicle"
    is_vehicle <- grepl("^veh$|^vehicle$", non_numeric_values)
    order_values[non_numeric_indices[is_vehicle]] <- min_num - 2  # Highest priority

    # Assign next priority to other non-numeric values
    order_values[non_numeric_indices[!is_vehicle]] <- min_num - 1
  }

  # Create a data frame for sorting
  df <- data.frame(
    original_value = column,
    order_value = rev(order_values),
    stringsAsFactors = FALSE
  )

  # Remove duplicates while preserving order
  df_unique <- df[!duplicated(df$original_value), ]

  # Sort the data frame by order_value
  df_ordered <- df_unique[order(df_unique$order_value), ]

  # Create an ordered factor with levels in the sorted order
  ordered_factor <- factor(
    column,
    levels = df_ordered$original_value,
    ordered = TRUE
  )

  return(ordered_factor)
}


