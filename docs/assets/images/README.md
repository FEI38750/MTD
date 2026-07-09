# MTD Explorer — website logo pack

This pack was generated directly from the approved logo image.

## Recommended files

### Homepage / hero
`mtd-explorer-logo.svg`

### MkDocs header
`mtd-explorer-mark.svg`

### Favicon
`favicon.ico`

### Apple icon
`apple-touch-icon.png`

### Social preview
`mtd-explorer-social-preview-1200x630.png`

## Suggested repository layout

Copy the pack into:

```text
docs/assets/images/
```

## MkDocs Material

```yaml
theme:
  name: material
  logo: assets/images/mtd-explorer-mark.svg
  favicon: assets/images/favicon.ico

extra:
  social:
    - icon: fontawesome/brands/github
      link: https://github.com/patrick-douglas/MTD
```

Use the full logo on the homepage:

```html
<img
  src="assets/images/mtd-explorer-logo.svg"
  alt="MTD Explorer logo"
  class="mtd-main-logo"
>
```

## Important note about the SVG files

The approved logo originated as a raster image. To preserve it exactly,
the SVG files in this pack contain the approved PNG embedded inside an SVG
container. They display identically and are convenient for the website, but
they are not a manually redrawn path-based vector illustration.

The 1254 × 1254 PNG master is more than sufficient for the homepage,
navigation header, retina screens, favicons, and normal browser zoom.
