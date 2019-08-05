using System;
using System.ComponentModel;
using System.Windows.Data;
using System.Windows.Media;
using QuickGraph;

namespace GraphSharp
{
    
    public class Vertex
    {
        public Object Name { get; set; }
        public Color Color { get; set; }
        private Color defaultColor = Colors.White;

        public Vertex(Object name, Color color)
            : base()
        {
            this.Name = name;
            this.Color = color;
        }

        public Vertex(Object name)
            : base()
        {
            this.Name = name;
            this.Color = defaultColor;
        }

        /// <summary>
        /// EXTENDED VERTEX METHODS
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>

        public bool Equals(Vertex other)
        {
            return other.Name == this.Name;
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((Vertex)obj);
        }

        public static bool operator==(Vertex v1, Vertex v2)
        {
            // If both are null, or both are same instance, return true.
            if (ReferenceEquals(v1, v2))
            {
                return true;
            }

            // If one is null, but not both, return false.
            if (((object)v1 == null) || ((object)v2 == null))
            {
                return false;
            }

            // Return true if the fields match:
            return v1.Name == v2.Name;
        }

        public static bool operator !=(Vertex v1, Vertex v2)
        {
            return !(v1 == v2);
        }

        public override int GetHashCode()
        {
            return (Name != null ? Name.GetHashCode() : 0);
        }

        public override string ToString()
        {
            return this.Name.ToString();
        }
    }
}