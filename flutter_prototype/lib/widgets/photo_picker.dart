import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import 'package:image_picker/image_picker.dart';
import 'dart:io';
import '../services/payment_service.dart';
import '../services/analysis_service.dart';

class PhotoPicker extends StatefulWidget {
  @override
  _PhotoPickerState createState() => _PhotoPickerState();
}

class _PhotoPickerState extends State<PhotoPicker> {
  final _urlController = TextEditingController();
  final _imagePicker = ImagePicker();
  
  @override
  Widget build(BuildContext context) {
    return Consumer2<PaymentService, AnalysisService>(
      builder: (context, paymentService, analysisService, child) {
        return Container(
          padding: EdgeInsets.all(40),
          decoration: BoxDecoration(
            color: Colors.white.withOpacity(0.08), // Interactive element background
            borderRadius: BorderRadius.circular(8),
          ),
          child: Column(
            mainAxisAlignment: MainAxisAlignment.center,
            children: [
              // Upload button
              _buildUploadButton(paymentService, analysisService),
              
              SizedBox(height: 20),
              
              // URL input section
              _buildUrlSection(paymentService, analysisService),
            ],
          ),
        );
      },
    );
  }
  
  Widget _buildUploadButton(PaymentService paymentService, AnalysisService analysisService) {
    return GestureDetector(
      onTap: () => _pickImage(paymentService, analysisService),
      child: Container(
        padding: EdgeInsets.all(40),
        decoration: BoxDecoration(
          border: Border.all(
            color: Colors.white.withOpacity(0.7),
            style: BorderStyle.solid,
            width: 1, // Minimal border per ui.mdc
          ),
          borderRadius: BorderRadius.circular(8),
        ),
        child: Column(
          children: [
            Icon(
              Icons.cloud_upload,
              size: 24,
              color: Colors.white.withOpacity(0.7),
            ),
            SizedBox(height: 12),
            Text(
              'Upload',
              style: TextStyle(
                color: Colors.white.withOpacity(0.7),
                fontSize: 16,
                fontWeight: FontWeight.w500,
              ),
            ),
          ],
        ),
      ),
    );
  }
  
  Widget _buildUrlSection(PaymentService paymentService, AnalysisService analysisService) {
    return Row(
      children: [
        Expanded(
          child: TextField(
            controller: _urlController,
            decoration: InputDecoration(
              hintText: 'Image URL...',
              isDense: true,
              contentPadding: EdgeInsets.symmetric(horizontal: 16, vertical: 12),
            ),
            style: TextStyle(color: Colors.white),
            onSubmitted: (_) => _analyzeUrl(paymentService, analysisService),
          ),
        ),
        SizedBox(width: 12),
        GestureDetector(
          onTap: () => _analyzeUrl(paymentService, analysisService),
          child: Container(
            padding: EdgeInsets.all(12),
            decoration: BoxDecoration(
              color: Colors.white.withOpacity(0.08),
              borderRadius: BorderRadius.circular(8),
            ),
            child: Icon(
              Icons.search,
              size: 16,
              color: Colors.white,
            ),
          ),
        ),
      ],
    );
  }
  
  Future<void> _pickImage(PaymentService paymentService, AnalysisService analysisService) async {
    // Check payment method
    if (!paymentService.canAnalyze()) {
      paymentService.showPayment();
      return;
    }
    
    try {
      final XFile? image = await _imagePicker.pickImage(
        source: ImageSource.gallery,
        maxWidth: 1920,
        maxHeight: 1080,
        imageQuality: 85,
      );
      
      if (image != null) {
        final imageFile = File(image.path);
        await analysisService.analyzeImage(imageFile);
        await paymentService.incrementUsage();
        
        ScaffoldMessenger.of(context).showSnackBar(
          SnackBar(content: Text('Image analysis complete')),
        );
      }
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Image pickup failed: $error')),
      );
    }
  }
  
  Future<void> _analyzeUrl(PaymentService paymentService, AnalysisService analysisService) async {
    final url = _urlController.text.trim();
    if (url.isEmpty) return;
    
    // Check payment method
    if (!paymentService.canAnalyze()) {
      paymentService.showPayment();
      return;
    }
    
    try {
      _urlController.clear();
      await analysisService.analyzeImageUrl(url);
      await paymentService.incrementUsage();
      
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('URL analysis complete')),
      );
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('URL analysis failed: $error')),
      );
    }
  }
  
  @override
  void dispose() {
    _urlController.dispose();
    super.dispose();
  }
} 