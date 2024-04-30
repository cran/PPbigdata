##Function to calculate CM index for nugget projection
CMNugg<- function(nuggproj,weight){
    1-HoleNugg(nuggproj,weight)
}
