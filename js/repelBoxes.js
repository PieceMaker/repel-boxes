//----------------------------------------------------------------------------------------------------------------------
/// RepelBoxes
///
/// @module
//----------------------------------------------------------------------------------------------------------------------

import normal from 'distributions-normal-random';
import _ from 'lodash';

class RepelBoxes {
    constructor() {
    } // end constructor

    toGeoJSONFeature(geom)
    {
        return {
            "type": "Feature",
            geometry: geom
        };
    } // end toGeoJSONFeature

    toGeoJSONPoint(x, y)
    {
        var point = {
            "type": "Point",
            "coordinates": [x, y]
        };
        return point;
    } // end toGeoJSONPoint

    addGeoJSONPoints(p1, p2)
    {
        var x = p1.coordinates[0] + p2.coordinates[0];
        var y = p1.coordinates[1] + p2.coordinates[1];
        return this.toGeoJSONFeature(
            this.toGeoJSONPoint(x, y)
        );
    } // end addGeoJSONPoints

    bboxToPolyCoords(bbox)
    {
        // bbox is of form [xmin, ymin, xmax, ymax]
        // polygon is of form [
        //     [xmin, ymin],
        //     [xmax, ymin],
        //     [xmax, ymax],
        //     [xmin, ymax]
        // ]
        return [
            [bbox[0], bbox[1]],
            [bbox[2], bbox[1]],
            [bbox[2], bbox[3]],
            [bbox[0], bbox[3]]
        ];
    } // end bboxToPolyCoords

    // Euclidean distance between two numeric vectors
    // p1    A 2D GeoJSON feature with point geometry
    // p2    A 2D GeoJSON feature with point geometry
    // Probably exists in turf.js
    // Returns a numeric
    euclid(p1, p2)
    {
        var retval = Math.pow(p1.geometry.coordinates[0] - p2.geometry.coordinates[0], 2);
        retval += Math.pow(p1.geometry.coordinates[1] - p2.geometry.coordinates[1], 2);
        return retval;
    } // end euclid

    // Find the intersections between a line and a rectangle
    // p1    A 2D GeoJSON point
    // p2    A 2D GeoJSON point
    // bbox  A 2D GeoJSON bounding box
    // Returns a point
    instersectLineRectangle(p1, p2, bbox)
    {
        var slope = (p2.coordinates[1] - p1.coordinates[1]) / (p2.coordinates[0] - p1.coordinates[0]);
        var intercept = p2.coordinates[1] - p2.coordinates * slope;
        var x, y;
        var retval = {
            "type": "GeometryCollection",
            "geometries": [
                {}, //Intersection 1
                {}, //Intersection 2
                {}, //Intersection 3
                {}  //Intersection 4
            ]
        };

        x = bbox[0];
        y = slope * x + intercept;
        if (bbox[1] <= y && y <= bbox[3]) {
            //TODO: create utility to create GeoJSON that is called like 'retval.geometries[0] = toGeoJSONPoint(x, y)'
            retval.geometries[0] = this.toGeoJSONPoint(x, y);
        }

        x = bbox[2];
        y = slope * x + intercept;
        if (bbox[1] <= y && y <= bbox[3]) {
            retval.geometries[1] = this.toGeoJSONPoint(x, y);
        }

        y = bbox[1];
        x = (y - intercept) / slope;
        if (bbox[0] <= y && y <= bbox[2]) {
            retval.geometries[2] = this.toGeoJSONPoint(x, y);
        }

        y = bbox[3];
        x = (y - intercept) / slope;
        if (bbox[0] <= y && y <= bbox[2]) {
            retval.geometries[3] = this.toGeoJSONPoint(x, y);
        }

        // Finding the nearest intersection on the rectangle so this function will only ever return a single
        // GeoJSON point
        var imin;
        var d = this.euclid(retval.geometries[0], p1);
        var dmin = d;
        for(var i = 1; i < 4; i++) {
            d = this.euclid(retval.geometries(i), p1);
            if(d < dmin) {
                dmin = d;
                imin = i;
            }
        }

        // Return nearest intersection point, if it exists
        // TODO: Return the point as GeoJSON
        return retval.geometries[i];
    } // end instersectLineRectangle

    // Move a box into the area specified by x limits and y limits
    // box    A 2D GeoJSON feature with polygon geometry and corresponding bbox key on geometry
    // bound  A 2D GeoJSON feature with polygon geometry and corresponding bbox key on geometry;
    //          representing the target bbox to move the first bbox to
    putWithinBounds(box, bound)
    {
        var d;
        if(box.bbox[0] < bound.bbox[0]) {
            d = Math.abs(box.bbox[0] - bound.bbox[0]);
            box.bbox[0] += d;
            box.bbox[2] += d;
        } else if(box.bbox[2] > bound.bbox[1]) {
            d = Math.abs(box.bbox[2] - bound.bbox[1]);
            box.bbox[0] -= d;
            box.bbox[2] -= d;
        }
        if(box.bbox[1] < bound.bbox[0]) {
            d = Math.abs(box.bbox[1] - bound.bbox[0]);
            box.bbox[1] += d;
            box.bbox[3] += d;
        } else if(box.bbox[3] > bound.bbox[1]) {
            d = Math.abs(box.bbox[3] - bound.bbox[1]);
            box.bbox[1] -= d;
            box.bbox[3] -= d;
        }

        box.geometry.coordinates = this.bboxToPolyCoords(box.bbox);

        return box;
    } // end putWithinBounds

    // Get the coordinates of the center of a box
    // bbox    The bbox key from a 2D GeoJSON feature
    // Returns a 2D GeoJSON feature with point geometry
    centroid(bbox)
    {
        var x = (bbox[0] + bbox[2]) / 2;
        var y = (bbox[1] + bbox[3]) / 2;
        return this.toGeoJSONFeature(
            this.toGeoJSONPoint(x, y)
        );
    } // end centroid

    // Test if a box overlaps another box
    // a  The bbox key from a 2D GeoJSON feature
    // b  The bbox key from a 2D GeoJSON feature
    // Returns a boolean
    overlaps(a, b)
    {
        return (b[0] <= a[2] &&
                b[1] <= a[3] &&
                b[2] >= a[0] &&
                b[3] >= a[1]);
    } // end overlaps

    // Test if a point is within the boundaries of a box
    // p    The coordinates key from a 2D GeoJSON feature
    // box  The bbox key from a 2D GeoJSON feature
    // Returns a boolean
    pointWithinBox(p, box)
    {
        return (p[0] >= box[0] &&
                p[0] <= box[2] &&
                p[1] >= box[1] &&
                p[1] <= box[3]);
    } // end pointWithinBox

    // Compute the repulsion force upon point a from point b.
    //
    // The force decays with the squared distance between the points,
    // similar to the force of repulsion between magnets.
    //
    // a      A 2D GeoJSON feature with point geometry
    // b      A 2D GeoJSON feature with point geometry
    // force  Magnitude of the force (defaults to 1e-6)
    // Returns a 2D GeoJSON feature with point geometry
    repelForce(a, b, force = 0.000001)
    {
        // Note that this is using a random value. Does this mean the label positions will converge
        // to different locations each time this is run?
        a.geometry.coordinates[0] += normal(1, {'mu': 0, 'sigma': force});
        a.geometry.coordinates[1] += normal(1, {'mu': 0, 'sigma': force});
        // Constrain the minimum distance to be at least 0.01
        var d = Math.max(this.euclid(a, b), 0.01);
        // Compute a unit vector in the direction of the force
        var v = this.toGeoJSONFeature(
            this.toGeoJSONPoint(
                (a.geometry.coordinates[0] - b.geometry.coordinates[0]) / d,
                (a.geometry.coordinates[1] - b.geometry.coordinates[1]) / d
            )
        );
        // Divide the force by the squared distance
        return this.toGeoJSONFeature(
            this.toGeoJSONPoint(
                v.geometry.coordinates[0] * force / Math.pow(d, 2),
                v.geometry.coordinates[1] * force / Math.pow(d, 2)
            )
        );
    } // end repelForce

    // Compute the spring force upon point a from point b
    //
    // The force increases with the distance between the points, similar
    // to Hooke's law for springs
    //
    // a      A 2D GeoJSON feature with point geometry
    // b      A 2D GeoJSON feature with point geometry
    // force  Magnitude of the force (defaults to 1e-6)
    // Returns a 2D GeoJSON feature with point geometry
    springForce(a, b, force = 0.000001)
    {
        var d = this.euclid(a, b);
        d = ((d < 0.01) ? 0 : d);
        // Compute a unit vector in the direction of the force
        var v = this.toGeoJSONFeature(
            this.toGeoJSONPoint(
                (a.geometry.coordinates[0] - b.geometry.coordinates[0]) / d,
                (a.geometry.coordinates[1] - b.geometry.coordinates[1]) / d
            )
        );
        return this.toGeoJSONFeature(
            this.toGeoJSONPoint(
                v.geometry.coordinates[0] * force * d,
                v.geometry.coordinates[1] * force * d
            )
        );
    } // end springForce

    // Adjust the layout of a list of potentially overlapping boxes
    //
    // centroids  A GeoJSON FeatureCollection where each feature has a bbox key representing the label bbox
    //              and each geometry is a feature point of the a label's centroid.
    // xlim       A numeric array representing the limits on the x axis like [xmin, xmax]
    // ylim       A numeric array representing the limits on the y axis like [ymin, ymax]
    // force      Magnitude of the force (defaults to 1e-6)
    // maxiter    Maximum number of iterations to try to resolve overlaps (defaults to 2000)
    repelBoxes(centroids, xlim, ylim, force = 0.000001, maxiter = 2000)
    {
        var i,j;
        var n = centroids.features.length();
        var iter = 0;
        var anyOverlaps = true;

        var ratios = [];
        var centroid;
        var originalCentroids = centroids;
        for(i = 0; i < n; i++) {
            centroid = centroids.features[i];
            // height over width
            ratios[i] = (centroids.bbox[3] - centroids.bbox[1]) / (centroids.bbox[2] - centroids.bbox[0]);
        }

        var totalForce = 1;
        var f;
        var ci;
        var cj;
        var fromBox, toBox;

        while(anyOverlaps && iter < maxiter && totalForce > 0) {
            if(iter == 0) {
                // Initialize
                totalForce = 0;
            }
            iter += 1;
            anyOverlaps = false;

            for(i = 0; i < n; i++) {
                f = this.toGeoJSONFeature(
                    this.toGeoJSONPoint(0, 0)
                );
                fromBox = centroids.features[i];
                ci = centroids.features[i];

                // Do we really want to compare all points to each other twice?
                // Maybe this could be changed to var j = i+1?
                for(j = 0; j < n; j++) {
                    toBox = centroids.features[j];
                    cj = centroids.features[j];

                    if(i == j) {
                        // Repel the box from its own centroid.
                        if(this.pointWithinBox(originalCentroids.features[i].coordinates, fromBox.bbox)) {
                            anyOverlaps = true;
                            f = this.addGeoJSONPoints(
                                f,
                                this.repelForce(ci, originalCentroids.features[i], force)
                            );
                        }
                    } else {
                        // Repel the box from overlapping boxes.
                        if(this.overlaps(fromBox, toBox)) {
                            anyOverlaps = true;
                            f = this.addGeoJSONPoints(
                                f,
                                this.repelForce(ci, cj, force)
                            );
                        }
                        // Repel the box from overlapping centroids.
                        if(this.pointWithinBox(originalCentroids.features[j].coordinates, fromBox.bbox)) {
                            anyOverlaps = true;
                            f = this.addGeoJSONPoints(
                                f,
                                this.repelForce(ci, originalCentroids.features[j], force)
                            );
                        }
                    }
                }

                // Pull toward the label's point
                if(!anyOverlaps) {
                    f = this.addGeoJSONPoints(
                        f,
                        this.springForce(originalCentroids.features[i], ci, force)
                    );
                }

                // Scale the x force by the ratio of height/width
                f.geometry.coordinates[0] = f.geometry.coordinates[0] * ratios[i];

                // Dampen the forces
                f.geometry.coordinates = f.geometry.coordinates * (1 - 0.001);

                totalForce += Math.abs(f.geometry.coordinates[0]) + Math.abs(f.geometry.coordinates[1]);

                centroids.features[i].bbox = [
                    fromBox.bbox[0] + f.geometry.coordinates[0],
                    fromBox.bbox[1] + f.geometry.coordinates[1],
                    fromBox.bbox[2] + f.geometry.coordinates[0],
                    fromBox.bbox[3] + f.geometry.coordinates[1]
                ];
                centroids.features[i] = this.putWithinBounds(centroids.features[i], xlim, ylim);
            }
        }

        for(i = 0; i < n; i++) {
            centroid = centroids.features[i];
            centroids.features[i] = this.toGeoJSONPoint(
                (centroid.bbox[0] + centroid.bbox[2]) / 2,
                (centroid.bbox[1] + centroid.bbox[3]) / 2
            );
        }

        return centroids;
    } // end repelBoxes;

} // end RepelBoxes

//----------------------------------------------------------------------------------------------------------------------

export default RepelBoxes;

//----------------------------------------------------------------------------------------------------------------------