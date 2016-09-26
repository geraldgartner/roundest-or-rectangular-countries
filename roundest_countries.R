library(sp)
library(rgeos)
library(maptools)
library(rgdal)

tau = 2*pi

world_shp = readShapePoly("~/Google Drive/dStd.at/2016-10-02-roundest_rectangular_countries/roundest-or-rectangular-countries/ne_10m_admin_0_sovereignty.shp")
world_shp = world_shp[which(world_shp$TYPE %in% c("Country", "Sovereign country", "Disputed")), ]

make_simple_poly = function(params, area) {
  # Function to return a rectangle as a SpatialPolygons object.
  # params is a 4-element vector.
  
  centroid = c(params[1], params[2])
  angle = params[3]
  aspect = params[4]
  
  lon_side_length = sqrt(area * aspect)
  lat_side_length = area / lon_side_length
  
  vertices = matrix(0, nrow=5, ncol=2)
  vertices[1, ] = c(-lon_side_length/2, -lat_side_length/2)
  vertices[2, ] = c(-lon_side_length/2, +lat_side_length/2)
  vertices[3, ] = c(lon_side_length/2, lat_side_length/2)
  vertices[4, ] = c(lon_side_length/2, -lat_side_length/2)
  vertices[5, ] = vertices[1, ]
  
  # Rotate:
  rot_mat = matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)),
                   nrow=2, ncol=2)
  
  vertices = t(rot_mat %*% t(vertices))
  
  # Centre:
  vertices[ , 1] = vertices[ , 1] + centroid[1]
  vertices[ , 2] = vertices[ , 2] + centroid[2]
  
  # Convert to SpatialPolygons:
  polygon_object = SpatialPolygons(list(Polygons(list(Polygon(vertices)), 0)))
  return(polygon_object)
}

num_features = length(world_shp$NAME)
overlaps = numeric(num_features)

# How many iterations before terminating/reducing the increment?
max_iter = 100

# Indices of the countries getting a second, third, ... run:
second_pass_i = c(26,28,31,41,44,45,53,55,59,63,66,72,80,93,96,102,116,118,123,124,135,139,140,
                  141,142,146,147,155,159,161,165,166,174,178,179,184,186,195,197,198,204,205,206)

third_pass_i = c(41,45,96,118,141,142,166)

fourth_pass_i = c(141,166)

fifth_pass_i = 141

second_pass = FALSE
third_pass = FALSE
fourth_pass = FALSE
fifth_pass = FALSE

if (fifth_pass) {
  fourth_pass = FALSE
  third_pass = FALSE
  second_pass = FALSE
}

if (fourth_pass) {
  third_pass = FALSE
  second_pass = FALSE
}

if (third_pass) {
  second_pass = FALSE
}

# Work in small batches to avoid memory problems:
for (j in 151:187) {
  if (fifth_pass) {
    i = fifth_pass_i[j]
  } else if (fourth_pass) {
    i = fourth_pass_i[j]
  } else if (third_pass) {
    i = third_pass_i[j]
  } else if (second_pass) {
    i = second_pass_i[j]
  } else {
    i = j
  }
  
  print(sprintf("Starting country %d: %s", i, world_shp$NAME[i]))
  iteration_printouts = character()
  
  # The country under study:
  this_shp = world_shp[i, ]
  
  # Get some initial parameters for the rectangle:
  centroid = gCentroid(this_shp)@coords
  bounding_box = gEnvelope(this_shp)
  bounding_box = bounding_box@bbox
  area = gArea(this_shp)
  aspect = (bounding_box[1, 2] - bounding_box[1, 1]) / (bounding_box[2, 2] - bounding_box[2, 1])
  
  angle = 0
  
  params = c(centroid[1],
             centroid[2],
             angle,
             aspect)
  
  
  # increments:
  inc = c((bounding_box[1, 2] - bounding_box[1, 1]) / 100,
          (bounding_box[2, 2] - bounding_box[2, 1]) / 100,
          1 * tau / 360,
          0.01)
  
  # Manual over-rides -- bit of a mess!
  if (world_shp$NAME[i] == "Bahamas") {
    params[1] = -78
    params[2] = 24.7
  }
  
  if (world_shp$NAME[i] == "Cape Verde") {
    params[1] = -23.6
    params[2] = 15.1
  }
  
  if (world_shp$NAME[i] == "Chile") {
    params[4] = 0.3
    inc[1] = 0.1
    inc[2] = 0.4
  }
  
  if (world_shp$NAME[i] == "Fiji") {
    params[1] = 178.5
    params[2] = -17.5
    params[4] = 1
    inc[1] = 0.03
    inc[2] = 0.03
  }
  
  if (world_shp$NAME[i] == "France") {
    params[1] = 2
    params[2] = 47
    params[4] = 1.2
    inc[1] = 0.1
    inc[2] = 0.1
  }
  
  if (world_shp$NAME[i] == "Kiribati") {
    params[1] = -157.4
    params[2] = 1.83
    params[4] = 1
    inc[1] = 0.004
    inc[2] = 0.004
  }
  
  if (world_shp$NAME[i] == "Maldives") {
    params[1] = 73.54
    params[2] = 1.91
    params[4] = .2
    inc[1] = 0.01
    inc[2] = 0.01
  }
  
  if (world_shp$NAME[i] == "Marshall Is.") {
    params[1] = 171.3
    params[2] = 7.07
    params[4] = 5
    inc[1] = 0.01
    inc[2] = 0.01
  }
  
  if (world_shp$NAME[i] == "Micronesia") {
    params[1] = 158.22
    params[2] = 6.87
    params[4] = 1
    inc[1] = 0.002
    inc[2] = 0.002
  }
  
  if (world_shp$NAME[i] == "Netherlands") {
    params[1] = 5.4
    params[2] = 52.3
    params[4] = 1
    inc[1] = 0.04
    inc[2] = 0.03
  }
  
  if (world_shp$NAME[i] == "New Zealand") {
    params[1] = 173
    params[2] = -41
    params[4] = 0.3
    inc[1] = 0.1
    inc[2] = 0.1
  }
  
  if (world_shp$NAME[i] == "Norway") {
    params[1] = 10
    params[2] = 64
    params[4] = 1.2
    inc[1] = 0.2
    inc[2] = 0.25
  }
  
  if (world_shp$NAME[i] == "Palau") {
    params[1] = 134.57
    params[2] = 7.47
    params[4] = 0.3
    inc[1] = 0.003
    inc[2] = 0.005
  }
  
  if (world_shp$NAME[i] == "Portugal") {
    params[1] = -8
    params[2] = 39.5
    params[4] = 0.4
    inc[1] = 0.02
    inc[2] = 0.05
  }
  
  if (world_shp$NAME[i] == "Seychelles") {
    params[1] = 55.45
    params[2] = -4.65
    params[4] = 1
    inc[1] = 0.0016
    inc[2] = 0.0016
  }
  
  if (world_shp$NAME[i] == "Tonga") {
    params[1] = -175.22
    params[2] = -21.19
    params[4] = 2
    inc[1] = 0.003
    inc[2] = 0.003
  }
  
  if (world_shp$NAME[i] == "Tuvalu") {
    params[1] = 178.68
    params[2] = -7.48
    params[4] = 1
    inc[1] = 0.0003
    inc[2] = 0.0003
  }
  
  if (world_shp$NAME[i] == "United States") {
    params[4] = 2.5
    inc[1] = 0.1
    inc[2] = 0.1
  }
  
  if (second_pass) {
    if (i == 26) { params = c(-54.130700, -8.323300, -0.532300, 1.570200) }
    if (i == 28) { params = c(114.742200, 4.569600, 0.619600, 1.819900) }
    if (i == 31) { params = c(19.877400, 6.399800, 0.340300, 2.138400) }
    if (i == 41) { params = c(43.364400, -11.660000, -0.509600, 1.009200) }
    if (i == 44) { params = c(-78.678300, 21.539400, -0.445100, 4.147500) }
    if (i == 45) { params = c(-68.906900, 12.137400, -0.277500, 2.334700) }
    if (i == 53) { params = c(-70.437700, 19.011200, -0.066300, 1.742900) }
    if (i == 55) { params = c(-78.419800, -1.392900, 0.710300, 1.413200) }
    if (i == 59) { params = c(25.675600, 58.722700, -0.115200, 2.218500) }
    if (i == 63) { params = c(2.774400, 46.455100, 0.089000, 1.350700) }
    if (i == 66) { params = c(-1.820000, 53.501700, 0.061100, 1.036600) }
    if (i == 72) { params = c(-14.983600, 12.051400, 0.136100, 2.193200) }
    if (i == 80) { params = c(-86.693800, 14.736700, 0.075000, 1.855300) }
    if (i == 93) { params = c(-77.196400, 18.148700, -0.132600, 2.994500) }
    if (i == 96) { params = c(137.698900, 36.251400, 0.574200, 2.478100) }
    if (i == 102) { params = c(-62.755900, 17.341400, -0.288000, 2.069500) }
    if (i == 116) { params = c(24.683900, 56.960200, -0.075000, 4.111100) }
    if (i == 118) { params = c(-5.823400, 32.177200, 0.151800, 1.362700) }
    if (i == 123) { params = c(-104.180200, 24.768500, -0.609300, 2.327900) }
    if (i == 124) { params = c(171.236000, 7.090000, -0.069800, 5.302000) }
    if (i == 135) { params = c(114.524600, 3.335800, 0.638800, 3.528600) }
    if (i == 139) { params = c(-85.073000, 12.905800, -0.146600, 0.820800) }
    if (i == 140) { params = c(5.672000, 52.384000, -0.548000, 0.782000) }
    if (i == 141) { params = c(10.180000, 62.425000, 0.144900, 1.353000) }
    if (i == 142) { params = c(83.517500, 28.341600, -0.427600, 3.099200) }
    if (i == 146) { params = c(68.142200, 29.013300, -0.254800, 1.038900) }
    if (i == 147) { params = c(-79.844500, 8.507300, 0.179800, 2.917300) }
    if (i == 155) { params = c(-58.339600, -23.318800, -0.745300, 2.093600) }
    if (i == 159) { params = c(104.253800, 61.455900, 0.022700, 7.751600) }
    if (i == 161) { params = c(43.904900, 24.468800, -0.740000, 1.743400) }
    if (i == 165) { params = c(103.817100, 1.359300, -0.090800, 2.285500) }
    if (i == 166) { params = c(160.090900, -8.613700, -0.289700, 3.038300) }
    if (i == 174) { params = c(19.588600, 48.692000, 0.127400, 4.176500) }
    if (i == 178) { params = c(-63.069800, 18.041200, -0.256600, 3.457300) }
    if (i == 179) { params = c(55.452100, -4.680600, 0.001700, 0.657000) }
    if (i == 184) { params = c(71.266100, 38.329300, -0.041900, 2.811400) }
    if (i == 186) { params = c(126.038200, -8.759600, 0.298500, 3.504200) }
    if (i == 195) { params = c(31.043200, 49.021300, -0.172800, 3.357800) }
    if (i == 197) { params = c(-101.485700, 39.638100, -0.176300, 2.679000) }
    if (i == 198) { params = c(62.195400, 41.915100, -0.389200, 3.111800) }
    if (i == 204) { params = c(-172.186800, -13.739800, -0.420600, 3.380400) }
    if (i == 205) { params = c(47.662600, 15.946700, 0.317600, 2.418300) }
    if (i == 206) { params = c(25.022800, -29.274300, 0.469500, 1.874900) }
    inc[4] = params[4]*0.01
  }
  
  if (third_pass) {
    if (i == 41) { params = c(43.360500, -11.644700, -0.354300, 0.808400) }
    if (i == 45) { params = c(-68.980800, 12.198700, -0.659700, 4.905200) }
    if (i == 96) { params = c(135.774000, 35.186000, 0.462500, 4.406100) }
    if (i == 118) { params = c(-7.597400, 30.973200, 0.596900, 2.863000) }
    if (i == 141) { params = c(10.600000, 62.550000, 0.254900, 1.622200) }
    if (i == 142) { params = c(83.981400, 28.223500, -0.373500, 3.994900) }
    if (i == 166) { params = c(159.971000, -8.505600, -0.361300, 3.752300) }
    inc[4] = params[4]*0.01
  }
  
  if (fourth_pass) {
    if (i == 141) { params = c(11.020000, 62.575000, 0.357900, 1.945000) }
    if (i == 166) { params = c(159.944400, -8.357600, -0.487000, 5.020600) }
    inc[4] = params[4]*0.01
  }
  
  if (fifth_pass) {
    # Fifth pass:
    # if (i == 141) { params = c(11.560000, 63.050000, 0.621400, 2.437100) }
    # # Bit hack-ish; only got Norway left now.
    # inc[1] = inc[1] * 0.1
    # inc[2] = inc[2] * 0.1
    # inc[4] = params[4]*0.03
    
    # Sixth pass:
    if (i == 141) { params = c(12.0040, 63.5150, 0.6493, 3.0366) }
    inc[1] = inc[1] * 0.1
    inc[2] = inc[2] * 0.1
    inc[4] = params[4] * 0.01
    
  }
  
  num_params = length(params)
  
  # this_rect will be the current rectangle
  this_rect = make_simple_poly(params, area)
  this_int = gIntersection(this_rect, this_shp)
  
  # In the case of a null intersection, gIntersection() returns NULL,
  # so we need to check !is.null before continuing (can't go ahead
  # and use gArea on a NULL).
  if (!is.null(this_int)) {
    this_overlap = gArea(this_int) / area
    
    # Counter for the iterations:
    iter_ct = 0
    
    # fine_incs will be TRUE after they get divided by 10:
    fine_incs = FALSE
    
    keep_iterating = TRUE
    
    # Initialising a couple of vectors for use in checking
    # if the algorithm is approximately converged:
    prev_prev_params = numeric(4)
    prev_params = numeric(4)
    
    while (keep_iterating) {
      if ((iter_ct == max_iter) && (!fine_incs)) {
        # Time to reduce the increments.
        inc = inc / 10
      }
      
      prev_overlap = this_overlap
      prev_prev_params = prev_params
      prev_params = params
      
      this_inc = numeric(num_params)
      
      for (i_param in 1:num_params) {
        # Increment the i-th parameter, see what the overlap is like:
        test_params = params
        test_params[i_param] = test_params[i_param] + inc[i_param]
        
        test_rect = make_simple_poly(test_params, area)
        test_int = gIntersection(test_rect, this_shp)
        
        if (!is.null(test_int)) {
          test_overlap = gArea(gIntersection(test_int, this_shp)) / area
          this_inc[i_param] = inc[i_param] * 2*(1*(test_overlap > prev_overlap) - 0.5)
        }
      }
      
      # Update the parameters:
      params = params + this_inc
      
      prev_rect = this_rect
      this_rect = make_simple_poly(params, area)
      this_overlap = gArea(gIntersection(this_rect, this_shp)) / area
      
      iteration_printouts = c(iteration_printouts,
                              sprintf("Step %d, overlap = %.3f, params = %.4f, %.4f, %.4f, %.4f",
                                      iter_ct, this_overlap, params[1], params[2], params[3], params[4]))
      
      print(iteration_printouts[iter_ct + 1])
      
      if (identical(params, prev_prev_params)) {
        if (fine_incs) {
          # Convergence!
          keep_iterating = FALSE
        } else {
          # Approximately converged, hopefully.
          inc = inc / 10
          fine_incs = TRUE
        }
      }
      
      
      iter_ct = iter_ct + 1
      if (iter_ct > 2*max_iter) {
        if (abs(this_overlap - prev_overlap) < 0.0005) {
          # I let it keep running a bit if it looks like
          # it's nowhere near converged.
          keep_iterating = FALSE
        } else {
          if (iter_ct > 5*max_iter) {
            keep_iterating = FALSE
          }
        }
      }
    }
    
    out_file = sprintf("images_no_title/%03d.svg", i)
    out_file_iter = sprintf("iterations/%03d.txt", i)
    
    # For plotting, want to have both the rectangle
    # and the country within the axis bounds.
    bbox_rect = gEnvelope(this_rect)@bbox
    min_x = min(bbox_rect[1, 1], bounding_box[1, 1])
    max_x = max(bbox_rect[1, 2], bounding_box[1, 2])
    min_y = min(bbox_rect[2, 1], bounding_box[2, 1])
    max_y = max(bbox_rect[2, 2], bounding_box[2, 2])
    
    # 303x346 should give approximately a 200x200 plotting area.
    svg(filename=out_file, width=303, height=346)
    plot(this_shp, xlim = c(min_x, max_x), ylim = c(min_y, max_y))
    plot(this_rect, add=TRUE)
    dev.off()
    
    out_conn = file(out_file_iter)
    writeLines(iteration_printouts, out_conn)
    close(out_conn)
    
    # Append the results to file:
    output_line = sprintf("%s: %.3f", world_shp$NAME[i], this_overlap)
    write(output_line, file="output.txt", append=TRUE)
    
    print(sprintf("%s: %.3f", world_shp$NAME[i], this_overlap))
  } else {
    # To be returned to later....
    print("No overlap!")
    print(world_shp$NAME[i])
    write(as.character(world_shp$NAME[i]), file="skipped_no_overlap.txt", append=TRUE)
  }
}