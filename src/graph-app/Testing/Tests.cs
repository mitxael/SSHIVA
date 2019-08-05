using System;
using System.Collections.Generic;
using System.Configuration;
using System.Linq;
using System.Text;
using CHoleR;
using GraphSharp;
using Microsoft.VisualStudio.TestTools.UnitTesting; ///Ref: Microsoft.QualityTools.UnitTestingFramework
using QuickGraph;

namespace CHoleR.Testing
{
    [TestClass]
    class Tests
    {
        public BidirectionalWeightedGraph<object, WeightedEdge<object>> G;

        [TestMethod]
        public void Setup()
        {
            G = new BidirectionalWeightedGraph<object, WeightedEdge<object>>();

            Assert.IsNotNull(G);
        }

        [TestMethod]
        public void Test01()
        {
            int edgesN = G.EdgeCount;

            Assert.IsTrue(edgesN > 0, "The actualCount was not greater than five");
        }
    }
}
