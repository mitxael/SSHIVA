using System;
using System.ComponentModel;
using System.Diagnostics;
using System.Windows.Data;
using System.Windows.Media;
using GraphSharp;
using QuickGraph;

namespace GraphSharp
{

public class WeightedEdge<TVertex>: Edge<TVertex>
	{
		public TVertex Source { get; set; }
		public TVertex Target { get; set; }
		public Double Weight { get; set; }
		public Color Color { get; set; }
	    private Color defaultColor = Colors.White;
		
		public WeightedEdge(TVertex source, TVertex target, Double weight, Color color)
			: base(source, target)
		{
			this.Source = source;
			this.Target = target;
			this.Weight = weight;
		    this.Color = color;
		}

		public WeightedEdge(TVertex source, TVertex target)
			: base(source, target)
		{
			this.Source = source;
			this.Target = target;
			this.Weight = 0.0;
		    this.Color = defaultColor;
		}

	    /// <summary>
	    /// EXTENDED EDGE METHODS
	    /// </summary>
	    /// <param name="other"></param>
	    /// <returns></returns>

	    public bool Equals(WeightedEdge<object> other)
	    {
	        return other.ToString() == this.ToString();
	    }
	    public override bool Equals(object obj)
	    {
	        if (ReferenceEquals(null, obj)) return false;
	        if (ReferenceEquals(this, obj)) return true;
	        if (obj.GetType() != this.GetType()) return false;
	        return Equals((WeightedEdge<TVertex>)obj);
	    }

		public static bool operator ==(WeightedEdge<TVertex> e1, WeightedEdge<TVertex> e2)
		{
			// If both are null, or both are same instance, return true.
			if (ReferenceEquals(e1, e2))
			{
				return true;
			}

			// If one is null, but not both, return false.
			if (((object)e1 == null) || ((object)e2 == null))
			{
				return false;
			}

			// Return true if the fields match:
			return e1.Source.ToString() == e2.Source.ToString() && e1.Target.ToString() == e2.Target.ToString();
		}

		public static bool operator !=(WeightedEdge<TVertex> e1, WeightedEdge<TVertex> e2)
		{
			return !(e1 == e2);
		}

		public bool Equals(WeightedEdge<TVertex> other)
		{
			return this.Source.ToString() == other.Source.ToString() && this.Target.ToString() == other.Target.ToString();
		}

		public override int GetHashCode()
		{
			unchecked
			{
				int hashCode = (Source != null ? Source.GetHashCode() : 0);
				hashCode = (hashCode * 397) ^ (Target != null ? Target.GetHashCode() : 0);
				hashCode = (hashCode * 397) ^ (int)Weight;
				return hashCode;
			}
		}

		public override string ToString()
		{
			return string.Format("({0},{1})={2}", this.Source, this.Target,this.Weight);
		}
	}
}