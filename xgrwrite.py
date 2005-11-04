colormap = [
    ("#ffffff","white"),     #  0
    ("#000000","black"),     #  1
    ("#ff0000","red"),       #  2
    ("#00ff00","green"),     #  3
    ("#0000ff","blue"),      #  4
    ("#ffff00","yellow"),    #  5
    ("#bc8f8f","brown"),     #  6
    ("#dcdcdc","grey"),      #  7
    ("#9400d3","violet"),    #  8
    ("#00ffff","cyan"),      #  9
    ("#ff00ff","magenta"),   # 10
    ("#ffa500","orange"),    # 11
    ("#7221bc","indigo"),    # 12
    ("#670748","maroon"),    # 13
    ("#40e0d0","turquoise"), # 14
    ("#008b00","green4"),    # 15
]

def getcolorid(col):
    if type(col) is str:
        return {
            "w": 0,
            "k": 1,
            "r": 2,
            "g": 3,
            "b": 4,
            "c": 9,
            "m": 10,
            "y": 5,
        }[str]
    else:
        rgbstr = "#%02x%02x%02x"%(int(col[0]*256),int(col[1]*256),int(col[2]*256))
        for i in range(len(colormap)):
            if colormap[i][0] == rgbstr:
                return i
        colormap.append((rgbstr,rgbstr))
        return len(colormap)-1
#     else:
#         raise "Unknown color format: "+repr(type(col))

def writecolormap(file):
    file.write("""\
    <!-- Color map -->
    <colormap>
""")
    for i in range(len(colormap)):
        file.write("""\
      <color-def id="%i" rgb="%s" name="%s"/>
"""%(i,colormap[i][0],colormap[i][1]))

    file.write("""\
    </colormap>
""")


# def colorstr(col):
#     colormap =
#     if type(col) is str:
#         return colormap[col]
#     else:
#         return "#%02x%02x%02x"%(int(col[0]*256),int(col[1]*256),int(col[2]*256))

#######################################


def writefig(file,fig):
    file.write("""\
<?xml version="1.0" standalone="no"?>
<!DOCTYPE grace SYSTEM "grace.dtd">
<grace version="59900" xmlns:grace="http://plasma-gate.weizmann.ac.il/Grace/">
  <!-- Description -->
  <description>
    <text/>
  </description>
  <!-- Page properties -->
  <page width="842" height="595" fill="yes" color-id="0"/>
  <!-- Data formats -->
  <data-formats>
    <dates reference="0" wrap="no" wrap-year="1950"/>
    <world format="%.8g"/>
  </data-formats>
""")

    writeframe(file,fig)

    file.write("""\
  <atext id="timestamp" active="no" offset="(0, 0)">
    <location x="0.03" y="0.03"/>
    <text-properties angle="0" justification="0">
      <face-spec font-id="0" color-id="1" char-size="1"/>
    </text-properties>
    <text><![CDATA[\${timestamp}]]></text>
    <text-frame type="0" offset="0">
      <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
      <fill-spec color-id="0" pattern-id="0"/>
    </text-frame>
    <pointer active="no">
      <arrow type="0" length="0" dl-ff="0" ll-ff="0"/>
    </pointer>
  </atext>
  <!-- Definitions -->
  <definitions>
""")

    writecolormap(file)

    file.write("""\
    <!-- Font map -->
    <fontmap>
      <font-def id="0" name="Times-Roman" fallback="Times-Roman"/>
      <font-def id="1" name="Times-Italic" fallback="Times-Italic"/>
      <font-def id="2" name="Times-Bold" fallback="Times-Bold"/>
      <font-def id="3" name="Times-BoldItalic" fallback="Times-BoldItalic"/>
      <font-def id="4" name="Helvetica" fallback="Helvetica"/>
      <font-def id="5" name="Helvetica-Oblique" fallback="Helvetica-Oblique"/>
      <font-def id="6" name="Helvetica-Bold" fallback="Helvetica-Bold"/>
      <font-def id="7" name="Helvetica-BoldOblique" fallback="Helvetica-BoldOblique"/>
      <font-def id="8" name="Courier" fallback="Courier"/>
      <font-def id="9" name="Courier-Oblique" fallback="Courier-Oblique"/>
      <font-def id="10" name="Courier-Bold" fallback="Courier-Bold"/>
      <font-def id="11" name="Courier-BoldOblique" fallback="Courier-BoldOblique"/>
      <font-def id="12" name="Symbol" fallback="Symbol"/>
      <font-def id="13" name="ZapfDingbats" fallback="ZapfDingbats"/>
    </fontmap>
    <!-- Size scales -->
    <scales font-size="0.028" line-width="0.0015"/>
  </definitions>
</grace>
""")



#######################################

def writeframe(file,fig):
    file.write("""\
  <frame id="frame" active="yes" type="0">
    <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
    <fill-spec color-id="0" pattern-id="0"/>
    <viewport xmin="0.15" xmax="1.25" ymin="0.15" ymax="0.85"/>
    <legend active="yes" length="0.04" vgap="0.01" hgap="0.01" single-symbol="no" invert="no">
      <face-spec font-id="0" color-id="1" char-size="1"/>
      <legframe anchor="(1, 1)" offset="(-0.05, -0.05)" justification="9">
        <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
        <fill-spec color-id="0" pattern-id="1"/>
      </legframe>
    </legend>
""")

    for a in fig.axes:
        writegraph(file,a)

    file.write("""\
    <atext id="title" active="no" offset="(0, 0.06)">
      <location x="0.5" y="1"/>
      <text-properties angle="0" justification="6">
        <face-spec font-id="0" color-id="1" char-size="1.5"/>
      </text-properties>
      <text><![CDATA[Title]]></text>
      <text-frame type="0" offset="0.005">
        <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
        <fill-spec color-id="0" pattern-id="0"/>
      </text-frame>
      <pointer active="no">
        <arrow type="0" length="1" dl-ff="1" ll-ff="0"/>
      </pointer>
    </atext>
    <atext id="subtitle" active="no" offset="(0, 0.02)">
      <location x="0.5" y="1"/>
      <text-properties angle="0" justification="6">
        <face-spec font-id="0" color-id="1" char-size="1"/>
      </text-properties>
      <text><![CDATA[Subtitle]]></text>
      <text-frame type="0" offset="0.005">
        <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
        <fill-spec color-id="0" pattern-id="0"/>
      </text-frame>
      <pointer active="no">
        <arrow type="0" length="1" dl-ff="1" ll-ff="0"/>
      </pointer>
    </atext>
  </frame>
""")


#######################################

def writegraph(file,axes):
    file.write("""\
    <graph id="graph" active="yes">
      <presentation-spec type="xy" stacked="no" bargap="0"/>
""" + """\
      <xscale min="%g" max="%g" type="Normal" invert="no"/>
"""%(axes.get_xlim()) + """\
      <yscale min="%g" max="%g" type="Normal" invert="no"/>
"""%(axes.get_ylim()) + """\
      <zscale norm="1"/>
      <locator type="1">
        <fixedpoint active="no" x="0" y="0"/>
        <xformat format="general" prec="6"/>
        <yformat format="general" prec="6"/>
      </locator>
      <axis id="x_axis" type="x" active="yes">
        <placement zero="no" offset="(0, 0)"/>
        <axisbar active="yes">
          <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
        </axisbar>
        <axislabel layout="parallel" offset="(0, 0.08)" side-placement="normal">
          <text-properties angle="0" justification="0">
            <face-spec font-id="0" color-id="1" char-size="1"/>
          </text-properties>
          <text/>
        </axislabel>
        <ticks major-step="0.2" minor-divisions="1" auto-ticking="6" rounded-position="yes">
          <userticks type="none"/>
          <tickmarks active="yes" side-placement="both">
            <major size="1" inout-placement="in" grid-lines="no">
              <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
            </major>
            <minor size="0.5" grid-lines="no" inout-placement="in">
              <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
            </minor>
          </tickmarks>
          <ticklabels active="yes" side-placement="normal" transform="" prepend="" append="" offset="auto" skip="0" stagger="0" start="auto" stop="auto">
            <text-properties angle="0" justification="0">
              <face-spec font-id="0" color-id="1" char-size="1"/>
            </text-properties>
            <format format="general" prec="5"/>
          </ticklabels>
        </ticks>
      </axis>
      <axis id="y_axis" type="y" active="yes">
        <placement zero="no" offset="(0, 0)"/>
        <axisbar active="yes">
          <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
        </axisbar>
        <axislabel layout="parallel" offset="(0, 0.08)" side-placement="normal">
          <text-properties angle="0" justification="0">
            <face-spec font-id="0" color-id="1" char-size="1"/>
          </text-properties>
          <text/>
        </axislabel>
        <ticks major-step="0.2" minor-divisions="1" auto-ticking="6" rounded-position="yes">
          <userticks type="none"/>
          <tickmarks active="yes" side-placement="both">
            <major size="1" inout-placement="in" grid-lines="no">
              <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
            </major>
            <minor size="0.5" grid-lines="no" inout-placement="in">
              <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
            </minor>
          </tickmarks>
          <ticklabels active="yes" side-placement="normal" transform="" prepend="" append="" offset="auto" skip="0" stagger="0" start="auto" stop="auto">
            <text-properties angle="0" justification="0">
              <face-spec font-id="0" color-id="1" char-size="1"/>
            </text-properties>
            <format format="general" prec="5"/>
          </ticklabels>
        </ticks>
      </axis>
""")

    for l in axes.lines:
        writeset(file,l)

    file.write("""\
    </graph>
""")

#######################################

def writeset(file,line):
    file.write("""\
      <set id="0x816f958" active="yes" type="xy" skip="0" skipmindist="0">
        <symbol type="0" size="1" char="65" font-id="0">
          <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
          <fill-spec color-id="1" pattern-id="0"/>
        </symbol>
        <line type="1" fill-type="0" fill-rule="winding" baseline-type="0" draw-baseline="no" draw-droplines="no">
""" + """\
          <line-spec color-id="%i" pattern-id="1" style-id="1" width="1"/>
"""%getcolorid(line.get_color()) + """\
          <fill-spec color-id="0" pattern-id="0"/>
        </line>
        <annotation active="no" type="2" prepend="" append="">
          <text-properties angle="0" justification="6">
            <face-spec font-id="0" color-id="1" char-size="1"/>
          </text-properties>
          <format format="general" prec="3"/>
        </annotation>
        <errorbar active="yes" side-placement="both">
          <barline size="1">
            <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
          </barline>
          <riserline arrow-clip="no" clip-length="0.1">
            <line-spec color-id="1" pattern-id="1" style-id="1" width="1"/>
          </riserline>
        </errorbar>
        <legend-entry>
          <text/>
        </legend-entry>
""")

    datapairs = zip(line.get_xdata(),line.get_ydata())
    file.write("""\
        <dataset cols="2" rows="%i" comment="data.dat">
"""%(len(datapairs),))

    for x,y in datapairs:
        file.write("""\
          <row X="%g" Y="%g"/>
"""%(x,y))
    file.write("""\
        </dataset>
""")

    file.write("""\
      </set>
""")

if __name__ == "__main__":
    import pylab

    pylab.plot([1,2,3],[1,3,2])

    file = open('tryout.xgr','w')
    writefig(file,pylab.gcf())
    file.close()

#    pylab.show()
