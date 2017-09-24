// NRT: Attenuating Sphere - to solve an approximation of sound uniformly emitting from a line

#pragma once

#include <math.h>
#include <limits>

namespace AttenuatingSphere
{
    // Helper structs
    struct Vector
    {
        double x, y, z;
    };

    struct LineSegment
    {
        Vector start, end;
    };

    struct Sphere
    {
        Vector center;
        double radius;
    };

    // Stores the result of finding the average attenuated direction.
    struct AverageAttenuatedDirection
    {
        Vector closest_point;
        Vector avg_direction;
        double closest_distance;
        double avg_spread;
    };

    // This is the main function to use.
    bool SolveAverageAttenuatedDirection(
        const Sphere& sphere,
        const LineSegment* lines,
        const uint32 num_lines,
        AverageAttenuatedDirection& result);

    // Helper functions
    double Dot(const Vector& lhs, const Vector& rhs)
    {
        return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
    }

    Vector operator+ (const Vector& lhs, const Vector& rhs)
    {
        Vector v = { lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
        return v;
    }

    Vector operator- (const Vector& lhs, const Vector& rhs)
    {
        Vector v = { lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
        return v;
    }

    Vector operator* (double c, const Vector& v)
    {
        Vector cv = { c*v.x, c*v.y, c*v.z };
        return cv;
    }

    Vector operator/ (const Vector& v, double c)
    {
        Vector cv = { v.x/c, v.y/c, v.z/c };
        return cv;
    }

    double LengthSquared(const Vector& v)
    {
        return Dot(v, v);
    }

    template<class T>
    const T& Max(const T& a, const T& b)
    {
        return (a < b) ? b : a;
    }

    Vector ClosestPointToLineSegment(
        const LineSegment& line_segement,
        const Vector& point)
    {    
        const Vector& direction = line_segement.end - line_segement.start;
        const double scale_along_segment = Dot(direction, point - line_segement.start);

        if (scale_along_segment <= 0.0)
        {
            return line_segement.start;
        }

        const double length_squared = LengthSquared(direction);

        if (scale_along_segment >= length_squared)
        {
            return line_segement.end;
        }

        return line_segement.start + (scale_along_segment / length_squared) * direction;
    }

    Vector ClosestPointToLine(
        const LineSegment& line_segement,
        const Vector& point)
    {
        const Vector direction = line_segement.end - line_segement.start;
        const double scale_along_segment = Dot(direction, point - line_segement.start);
        const double length_squared = LengthSquared(direction);
        return line_segement.start + (scale_along_segment / length_squared) * direction;
    }

    // returns a LineSegment with endpoints 
    // assumes LineSegment is valid with two distinct end points
    bool ClipLineWithSphere(
        const Sphere& sphere,
        LineSegment& input)
    {        
        const Vector closest_to_line = ClosestPointToLine(input, sphere.center);
        const double distance_to_line_sq = LengthSquared(closest_to_line - sphere.center);
        const double radius_sq = sphere.radius * sphere.radius;
        if (distance_to_line_sq + DBL_EPSILON >= radius_sq)
        {
            return false;
        }
        
        const Vector direction = input.end - input.start;
        const double segment_length_sq = LengthSquared(direction);
        const double half_chord_length_sq = radius_sq - distance_to_line_sq;
        const double length_to_intersect = sqrt(half_chord_length_sq / segment_length_sq);
        
        const Vector intersect = length_to_intersect * direction;

        if (LengthSquared(input.start - sphere.center) > radius_sq)
        {
            input.start = closest_to_line - intersect;
        }

        if (LengthSquared(input.end - sphere.center) > radius_sq)
        {
            input.end = closest_to_line + intersect;
        }
        
        return true;
    }

    // Analytical solution for the total average, weighted direction a line segment.
    Vector TotalAttenuatedDirection(
        const Sphere& sphere,
        const LineSegment& line,
        double& out_spread)
    {
        LineSegment clipped_line(line);
        if (ClipLineWithSphere(sphere, clipped_line))
        {
            const double inv_radius = 1.0 / sphere.radius;
            const Vector start = inv_radius * (clipped_line.start - sphere.center);
            const Vector direction = inv_radius * (clipped_line.end - clipped_line.start);
            
            // Solve the line integral: integrate ((i/sqrt(distance to point)) * W)dt over the parametric form vStart + t*vDirection.
            // These are integrals of the form sqrt(ax^2 + bx + c):
            // https://en.wikipedia.org/wiki/List_of_integrals_of_irrational_functions#Integrals_involving_R_.3D_.E2.88.9Aax2_.2B_bx_.2B_c
            // which a, b, and c can be expressed as dot products.
            const double dot_end = LengthSquared(direction); // a
            if (dot_end < DBL_EPSILON)
            {
                return{ 0.0, 0.0, 0.0 };
            }
            const double dot_start = LengthSquared(start); // b
            const double dot2 = 2 * Dot(start, direction); // c

            const double poly1 = sqrt(dot_end + dot2 + dot_start); // distance at 1
            const double poly0 = sqrt(dot_start); // distance at 0

            const double length = sqrt(dot_end);

            // integral of the inverse distance
            const double int_start1 = log(Max(DBL_EPSILON, 2 * length * poly1 + 2 * dot_end + dot2)) / length;
            const double int_start0 = log(Max(DBL_EPSILON, 2 * length * poly0 + dot2)) / length;

            // integral of the normalized vector length
            const double int_end1 = (poly1 / dot_end) - (dot2 * int_start1) / (2 * dot_end);
            const double int_end0 = (poly0 / dot_end) - (dot2 * int_start0) / (2 * dot_end);

            // Spread
            {
                // Arc length is the dot product of the normalized average direction and the vector to the point.
                // Taking the limit of the dot product to find the total arc length simplifies to the length of the total direction.
                // To find the average spread needs to be normalized by the total weight. This is the line integral of the form:
                // integrate (W)dt or integrate (1 - sqrt(distance to point))dt.
                const double discriminent = 4 * dot_end * dot_start - dot2 * dot2;
                const double total_weighting1 = 1 - ((2 * dot_end + dot2) * poly1) / (4 * dot_end) -
                    (discriminent * int_start1) / (8 * dot_end);
                const double total_weighting0 = -(dot2 * poly0) / (4 * dot_end) -
                    (discriminent * int_start0) / (8 * dot_end);

                out_spread = length * (total_weighting1 - total_weighting0);
            }

            // definite integrals
            const double definite_start = int_start1 - int_start0;
            const double definite_end = int_end1 - int_end0;

            const double x = start.x * definite_start + direction.x * definite_end -
                start.x - direction.x / 2.0;
            const double y = start.y * definite_start + direction.y * definite_end -
                start.y - direction.y / 2.0;
            const double z = start.z * definite_start + direction.z * definite_end -
                start.z - direction.z / 2.0;

            const Vector total_direction = { x, y, z };
            return length * total_direction;
        }
        return { 0.0, 0.0, 0.0 };
    }

    bool SolveAverageAttenuatedDirection(
        const Sphere& sphere,
        const LineSegment* lines,
        const uint32 num_lines,
        AverageAttenuatedDirection& result)
    {
        Vector total_direction = { 0.0, 0.0, 0.0 };
        double total_weighting = 0.0;

        double closest_distance_sq = DBL_MAX;
        bool within_range = false;
        const double radius_sq = sphere.radius * sphere.radius;

        for (const LineSegment* line = lines; line != lines + num_lines; ++line)
        {
            Vector closest_point = ClosestPointToLineSegment(*line, sphere.center);
            const double distance_sq = LengthSquared(closest_point - sphere.center);

            if (distance_sq < closest_distance_sq)
            {
                closest_distance_sq = distance_sq;
                result.closest_point = closest_point;
            }

            if (distance_sq + DBL_EPSILON < radius_sq)
            {
                double weighting = 0;
                total_direction = total_direction + TotalAttenuatedDirection(sphere, *line, weighting);
                total_weighting += weighting;
                within_range = true;
            }
        }

        double total_length = sqrt(LengthSquared(total_direction));
        result.closest_distance = sqrt(closest_distance_sq);
                        
        if (total_length > DBL_EPSILON)
        {
            result.avg_spread = 1.0 - total_length / total_weighting;
            result.avg_direction = total_direction / total_length;
        }
        // If length went to zero, the distribution of lines cancelled out, for example two parallel lines.
        else if (closest_distance_sq < radius_sq)
        {
            result.avg_spread = 1.0f;
            result.avg_direction = result.closest_point;
        }
        else
        {
            result.avg_spread = 0.0f;
            result.avg_direction = result.closest_point;
        }

        return within_range;
    }
}