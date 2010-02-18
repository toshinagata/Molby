#
#  transform.rb
#
#  Created by Toshi Nagata.
#  Copyright 2008 Toshi Nagata. All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

class Transform

  #  Transform matrix that translates cx to origin and rotates ax->x
  #  and ay->y. Ax and ay will be normalized, and if ax and ay are
  #  not perpendicular, (ax.cross(ay)).cross(ax).normalize is used
  #  instead of ay.
  def self.rotation_with_axis(ax, ay, cx = Vector3D[0,0,0])
    ax = ax.normalize
    az = ax.cross(ay).normalize
    ay = az.cross(ax).normalize
    Transform[ax, ay, az, [0,0,0]].inverse * Transform.translation(cx * (-1))
  end

  #  Decompose a rotation matrix to the rotation axis and angle.
  #  Returns a list [axis, angle]
  def to_rot
    xx = self[0,0]
    yy = self[1,1]
    zz = self[2,2]
    ww = xx + yy + zz + 1.0
    if (ww >= 1.0)
      w = Math.sqrt(0.25 * ww)
      ww = 0.25 / w;
      x = (self[1,2] - self[2,1]) * ww
      y = (self[2,0] - self[0,2]) * ww
      z = (self[0,1] - self[1,0]) * ww
    elsif (xx >= yy && xx >= zz)
      x = Math.sqrt(0.25 * (1.0 + xx - yy - zz))
      ww = 0.25 / x
      y = (self[1,0] + self[0,1]) * ww
      z = (self[2,0] + self[0,2]) * ww
      w = (self[1,2] - self[2,1]) * ww
    elsif (yy >= zz)
      y = Math.sqrt(0.25 * (1.0 + yy - xx - zz))
      ww = 0.25 / y
      x = (self[1,0] + self[0,1]) * ww
      z = (self[1,2] + self[2,1]) * ww
      w = (self[2,0] - self[0,2]) * ww
    else
      z = Math.sqrt(0.25 * (1.0 + zz - xx - yy))
      ww = 0.25 / z
      x = (self[2,0] + self[0,2]) * ww
      y = (self[2,1] + self[1,2]) * ww
      w = (self[0,1] - self[1,0]) * ww
    end
#    puts "#{x}, #{y}, #{z}, #{w}"
    ww = Math.sqrt(x*x + y*y + z*z)
    ang = Math.atan2(ww, w) * 2
    v = Vector3D[x/ww, y/ww, z/ww]
    if (ang > Math::PI)
      ang = 2*Math::PI - ang
      v = v * (-1)
    end
    return [v, ang * Rad2Deg]
  end

  #  Decompose a rotation matrix to spiral components
  #  Returns a list [axis, angle, center, pitch]
  def to_spiral
    r, ang = to_rot   #  axis and angle
    d = column(3)     #  translation vector
    #  Calculate psuedoinverse (Moore-Penrose inverse) matrix of (I-R)
    #  and solve d = (I-R)c + (d.r)*r for c
    r0 = Vector3D[1,0,0].cross(r)
    if r0.length < 1e-4
      r0 = Vector3D[0,1,0].cross(r)
    end
    r0 = r0.normalize
    r1 = r.cross(r0).normalize
    q0 = r0 - (self * r0 - d)
    k0 = q0.length
    q0 = q0.normalize
    q1 = r1 - (self * r1 - d)
    k1 = q1.length
    q1 = q1.normalize
    u = Transform[q0, q1, r, [0,0,0]]
    uinv = u.inverse
    udet = u.determinant
#    puts "u = #{u}, udet = #{udet}, uinv = #{uinv}"
    usv = u * Transform.diagonal(k0, k1, 0) * (Transform[r0, r1, r, [0,0,0]].inverse)
#    puts "usv = #{usv}, rot = #{rot}"
    mpinv = Transform[r0 * (1.0/k0), r1 * (1.0/k1), [0,0,0], [0,0,0]] * (Transform[q0, q1, r, [0,0,0]].inverse)
#    mm = ir - ir * mpinv * ir
#    mm2 = mpinv - mpinv * ir * mpinv
#    puts "mpinv = #{mpinv}, mm = #{mm}, mm2 = #{mm2}"
    pitch = d.dot(r)
    c = mpinv * (d - r * pitch)
#    c2 = (rot2 * c - c).cross(r)
#    puts "c = #{c}, c2 = #{c2}"
    return [r, ang, c, pitch]
  end
end
