# Project Principles & Requirements

## Core Philosophy
- **Keep UI stupid simple**: No extraneous lines or outlining
- **Focus on code rather than markdown**: Don't overgenerate, focus on important aspects and correctness
- **Git commit with message after acceptance**: Every change must be committed

## Molecular Visualization
- **Sphere representation**: Use spheres that are of Bohr radius (van der Waals radii)
- **NOT ball-and-stick**: Avoid traditional ball-and-stick molecular models
- **Scale factor**: 0.8 gives good visual balance for van der Waals radii
- **Background**: Black (#000000) to match app theme
- **3D Library**: Use $3Dmol.js with sphere representation

## UI Design Principles
- **Minimal visual elements**: No unnecessary borders, shadows, or decorative elements
- **Functional design**: Every element serves a purpose
- **Clean typography**: Enhanced readability without visual clutter
- **Mobile-first**: Consider mobile device constraints and touch events
- **Progressive enhancement**: Gracefully degrade if features aren't available

## Technical Requirements
- **Modern JavaScript**: ES6+ features, async/await, modern DOM APIs
- **Error handling**: Always include proper try-catch blocks for async operations
- **HTTPS requirement**: Camera API requires HTTPS on mobile devices
- **Response format**: Consistent JSON response format with `{ output: {...} }` structure

## Molecular Data Standards
- **SMILES format only**: Accept valid SMILES notation, reject molecular formulas like "CaCO3"
- **SDF generation**: Always validate SMILES before attempting 3D coordinate generation
- **File organization**: Store SDF files in dedicated `sdf_files/` directory
- **Batch processing**: Handle multiple molecules efficiently in arrays

## User Experience
- **Immediate feedback**: Provide visual feedback for user interactions (loading states)
- **Column layout**: Add new object columns to the right for multiple analysis results
- **Descriptive names**: Show molecule names above 3D viewers for better context
- **Individual controls**: Allow closing individual analysis columns

## AI Integration
- **OpenAI Vision API**: Handle connection errors gracefully
- **Fallback behavior**: Provide fallback when object identification fails
- **SMILES preference**: Prefer valid SMILES over molecular formulas
- **Dual analysis**: Use both full image and cropped region for enhanced analysis

## Development Workflow
- **Live reload**: Maintain nodemon configuration for efficient development
- **SSL certificates**: Use existing certificates for HTTPS development
- **Debug configuration**: Leverage VS Code's auto-attach debugging feature
- **No comments**: Avoid adding comments unless absolutely necessary for complex algorithms 